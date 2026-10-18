"""
chain_dynamics.py
=================
Parse and plot grain-chain position data written as whitespace-separated
``.dat`` files (one row per grain, ``n_grain`` rows per time step, columns
x y z), as produced by the LSDEM simulations / experiment tracking used in
YL_postprocess.py.

Key variables (all user-settable, nothing hard-coded to the example file):
    data_dir   : folder that holds the .dat files (the database lives elsewhere)
    filename   : file inside data_dir, e.g. "positions_406.dat"
    n_grain    : number of grains in the chain -- "auto" (default) detects it
                 from the data, or pass an int to set it manually
    drive_end  : which end is driven -- "first" (grain 0 in the file) or
                 "last" (grain n_grain-1 in the file)

Grain numbering convention used by every analysis/plot function:
    grain 0 = the DRIVEN grain, grain n_grain-1 = the far (free) end.
If drive_end == "last", the grain axis is flipped internally so this always
holds; ``ChainData.raw_index(g)`` maps back to the row order in the file.

Quick use (Python):
    from chain_dynamics import ChainData
    chain = ChainData.from_file(r"D:\\LSDEM_DataCenter\\...\\Output",
                                "positions_406.dat", n_grain="auto",  # or 19
                                drive_end="first", frame_rate=240.0)
    chain.plot_snapshot(at=300)                       # every grain at one time
    chain.plot_time_history(grains=[0, 5, 10], start=0, end=600)

Quick use (command line):
    python chain_dynamics.py --data-dir D:\\path\\to\\Output --file positions_406.dat \
        --n-grain 19 --drive-end first --frame-rate 240 --snapshot 300 \
        --window 0 600 --grains 0 5 10 --save-dir figs
"""
from __future__ import annotations

import argparse
import glob
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

COMPONENTS = {"x": 0, "y": 1, "z": 2}
_FIRST_ALIASES = {"first", "start", "0", "grain0", "head"}
_LAST_ALIASES = {"last", "end", "-1", "tail"}


# ----------------------------------------------------------------------------
# File-level helpers
# ----------------------------------------------------------------------------
def detect_n_grain(rows: np.ndarray, max_lag: int | None = None,
                   max_sample: int = 200_000, n_probe: int = 2000,
                   dip: float = 0.5, return_scores: bool = False):
    """Auto-detect the number of grains per time step from raw rows (n_rows, n_cols).

    Idea: row i and row i+n_grain are the SAME grain one output step apart, so
    their distance is small (motion per frame); rows at any other lag compare
    different grains (distance ~ grain spacing). Score each lag by the median
    distance ||row[i+lag] - row[i]||; the true n_grain is the SMALLEST lag
    k >= 2 that is a sharp dip relative to its neighbours k-1 and k+1
    (score[k] < dip * min(score[k-1], score[k+1])), confirmed by a second dip
    at 2k. Taking the smallest dip (not the global minimum) avoids picking a
    multiple of n_grain when the motion is periodic in time.

    Returns n_grain (and the computed score array, indexed by lag, if return_scores).
    Raises ValueError if no dip is found (e.g. only one time step in the file,
    or grains move more than one spacing per output step).
    """
    rows = np.asarray(rows, dtype=np.float64)[:max_sample]
    n_rows = rows.shape[0]
    if max_lag is None:
        max_lag = n_rows // 3          # need room for the 2k confirmation
    max_lag = min(max_lag, (n_rows - 2) // 2)
    if max_lag < 2:
        raise ValueError("not enough rows to detect n_grain; pass n_grain manually")

    # fixed set of probe rows so each lag costs O(n_probe), not O(n_rows)
    last_start = n_rows - (2 * max_lag + 2)
    probe = np.unique(np.linspace(0, last_start, min(n_probe, last_start + 1)).astype(int))
    base = rows[probe]
    scores = np.full(2 * max_lag + 2, np.nan)

    def score(lag):
        if np.isnan(scores[lag]):
            scores[lag] = np.median(np.linalg.norm(rows[probe + lag] - base, axis=1))
        return scores[lag]

    def is_dip(k):
        return score(k) < dip * min(score(k - 1), score(k + 1))

    for k in range(2, max_lag + 1):
        if is_dip(k) and is_dip(2 * k):
            computed = scores[: 2 * k + 2]
            return (k, computed) if return_scores else k
    raise ValueError(f"could not auto-detect n_grain (no periodic dip for lags 2..{max_lag}); "
                     "pass n_grain manually")


def load_positions(file_path: str, n_grain: int | str | None = None, n_cols: int = 3) -> np.ndarray:
    """Read a positions .dat file into an array of shape (n_steps, n_grain, n_cols).

    Rows are grains, every ``n_grain`` consecutive rows form one time step.
    n_grain : int to set it manually, or None / 'auto' to auto-detect it
              (see detect_n_grain). If given manually and the auto-detected
              value disagrees, a warning is issued but the manual value is used.
    An incomplete trailing step (e.g. a simulation that was killed mid-write)
    is dropped with a warning instead of raising.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(file_path)

    data = pd.read_csv(file_path, sep=r"\s+", header=None,
                       dtype=np.float64, comment="#").values
    if data.shape[1] != n_cols:
        raise ValueError(f"{file_path}: expected {n_cols} columns, found {data.shape[1]}")

    name = os.path.basename(file_path)
    if n_grain is None or str(n_grain).lower() == "auto":
        n_grain = detect_n_grain(data)
        print(f"{name}: auto-detected n_grain = {n_grain}")
    else:
        n_grain = int(n_grain)
        if n_grain <= 0:
            raise ValueError("n_grain must be a positive integer")
        try:
            detected = detect_n_grain(data)
            if detected != n_grain:
                warnings.warn(f"{name}: n_grain={n_grain} given, but data looks like "
                              f"{detected} grains per step. Using {n_grain}.")
        except ValueError:
            pass  # detection inconclusive; trust the user

    n_steps, remainder = divmod(data.shape[0], n_grain)
    if n_steps == 0:
        raise ValueError(f"{file_path}: fewer rows ({data.shape[0]}) than n_grain ({n_grain})")
    if remainder:
        warnings.warn(f"{name}: {data.shape[0]} rows is not a multiple of "
                      f"n_grain={n_grain}; dropping the last {remainder} rows. "
                      "Check n_grain if this is unexpected.")
    return data[: n_steps * n_grain].reshape(n_steps, n_grain, n_cols)


def find_data_files(data_dir: str, pattern: str = "positions_*.dat") -> list[str]:
    """List matching data files in ``data_dir`` (sorted)."""
    if not os.path.isdir(data_dir):
        raise NotADirectoryError(data_dir)
    return sorted(glob.glob(os.path.join(data_dir, pattern)))


def normalize_drive_end(drive_end) -> str:
    """Accept 'first'/'last', 0/-1, 'start'/'end', ... and return 'first' or 'last'."""
    key = str(drive_end).strip().lower()
    if key in _FIRST_ALIASES:
        return "first"
    if key in _LAST_ALIASES:
        return "last"
    raise ValueError(f"drive_end must be 'first' (grain 0) or 'last' (end grain), got {drive_end!r}")


# ----------------------------------------------------------------------------
# Main container
# ----------------------------------------------------------------------------
class ChainData:
    """Positions of a grain chain over time, indexed so grain 0 is the driven grain.

    Parameters
    ----------
    positions   : array (n_steps, n_grain, 3) in FILE order.
    drive_end   : 'first' or 'last' -- which file-order end is driven.
    frame_rate  : output samples per second (experiments, e.g. 240 Hz).
    dt_out      : seconds between output steps (simulations: dt*time_scale*nTout).
                  Give frame_rate OR dt_out; if neither, time is in output steps.
    length_unit : label used on axes (e.g. 'mm', 'voxel'); data is not rescaled.
    length_scale: multiply raw coordinates by this when building displacements
                  (e.g. 4e-4 to convert voxel -> m). Default 1.0 (no change).
    """

    def __init__(self, positions: np.ndarray, drive_end="first", frame_rate=None,
                 dt_out=None, length_unit="mm", length_scale=1.0, name=""):
        positions = np.asarray(positions, dtype=np.float64)
        if positions.ndim != 3 or positions.shape[2] != 3:
            raise ValueError("positions must have shape (n_steps, n_grain, 3)")
        if frame_rate is not None and dt_out is not None:
            raise ValueError("give frame_rate OR dt_out, not both")

        self.drive_end = normalize_drive_end(drive_end)
        self.raw_positions = positions * length_scale
        self.positions = (self.raw_positions if self.drive_end == "first"
                          else self.raw_positions[:, ::-1, :])
        self.length_unit = length_unit
        self.name = name
        if frame_rate is not None:
            self.dt = 1.0 / float(frame_rate)
        elif dt_out is not None:
            self.dt = float(dt_out)
        else:
            self.dt = None  # time axis in steps

    # ---- construction --------------------------------------------------------
    @classmethod
    def from_file(cls, data_dir: str, filename: str, n_grain="auto", drive_end="first", **kwargs):
        """Load data_dir/filename. n_grain: int (manual) or 'auto'/None (auto-detect)."""
        path = os.path.join(data_dir, filename)
        kwargs.setdefault("name", os.path.splitext(filename)[0])
        return cls(load_positions(path, n_grain), drive_end=drive_end, **kwargs)

    # ---- basic properties ------------------------------------------------------
    @property
    def n_steps(self) -> int:
        return self.positions.shape[0]

    @property
    def n_grain(self) -> int:
        return self.positions.shape[1]

    @property
    def times(self) -> np.ndarray:
        """Time of every output step (s), or step index if no time base was given."""
        steps = np.arange(self.n_steps, dtype=float)
        return steps * self.dt if self.dt is not None else steps

    @property
    def time_label(self) -> str:
        return "Time (s)" if self.dt is not None else "Output step"

    def raw_index(self, g: int) -> int:
        """Drive-relative grain index -> row index inside each step of the file."""
        g = int(g) % self.n_grain
        return g if self.drive_end == "first" else self.n_grain - 1 - g

    def summary(self) -> str:
        span = f"{self.times[-1]:.4g} s" if self.dt is not None else f"{self.n_steps} steps"
        rep = self.repeated_frames()
        return (f"{self.name or 'chain'}: {self.n_steps} steps x {self.n_grain} grains, "
                f"driven end = {self.drive_end} (file order), duration = {span}, "
                f"repeated (identical) frames = {len(rep)}")

    def repeated_frames(self) -> np.ndarray:
        """Steps whose positions are identical to the previous step (dropped/duplicated frames)."""
        same = np.all(self.positions[1:] == self.positions[:-1], axis=(1, 2))
        return np.nonzero(same)[0] + 1

    # ---- index / time resolution --------------------------------------------
    def to_step(self, value, unit="step") -> int:
        """Convert a time (unit='s') or step (unit='step') to a valid step index."""
        if unit == "step":
            step = int(value)
        elif unit == "s":
            if self.dt is None:
                raise ValueError("unit='s' needs frame_rate or dt_out")
            step = int(round(float(value) / self.dt))
        else:
            raise ValueError("unit must be 'step' or 's'")
        if step < 0:
            step += self.n_steps
        if not 0 <= step < self.n_steps:
            raise IndexError(f"step {step} outside 0..{self.n_steps - 1}")
        return step

    def window(self, start=None, end=None, unit="step") -> slice:
        """Slice of steps [start, end] inclusive; None = from beginning / to end."""
        s = 0 if start is None else self.to_step(start, unit)
        e = self.n_steps - 1 if end is None else self.to_step(end, unit)
        if e < s:
            raise ValueError(f"window end ({e}) is before start ({s})")
        return slice(s, e + 1)

    def _grains(self, grains) -> list[int]:
        if grains is None:
            return list(range(self.n_grain))
        if np.isscalar(grains):
            grains = [grains]
        out = []
        for g in grains:
            g = int(g)
            if not -self.n_grain <= g < self.n_grain:
                raise IndexError(f"grain {g} outside 0..{self.n_grain - 1}")
            out.append(g % self.n_grain)
        return out

    # ---- data products ------------------------------------------------------
    def coordinate(self, component="z") -> np.ndarray:
        """Absolute coordinate, shape (n_steps, n_grain)."""
        return self.positions[:, :, COMPONENTS[component]]

    def displacement(self, component="z", reference="initial") -> np.ndarray:
        """Displacement of every grain, shape (n_steps, n_grain).

        component : 'x', 'y', 'z', or 'mag' (|vector displacement|).
        reference : 'initial' -> relative to step 0
                    'mean'    -> relative to each grain's time mean (as in YL_postprocess)
                    int       -> relative to that step
                    'none'    -> absolute coordinate (no subtraction)
        """
        pos = self.positions
        if reference == "none":
            if component == "mag":
                raise ValueError("component='mag' needs a reference")
            return pos[:, :, COMPONENTS[component]].copy()
        if reference == "initial":
            ref = pos[0]
        elif reference == "mean":
            ref = pos.mean(axis=0)
        elif isinstance(reference, (int, np.integer)):
            ref = pos[self.to_step(reference)]
        else:
            raise ValueError("reference must be 'initial', 'mean', 'none' or a step index")
        disp = pos - ref[np.newaxis]
        if component == "mag":
            return np.linalg.norm(disp, axis=2)
        return disp[:, :, COMPONENTS[component]]

    def _ylabel(self, component, reference) -> str:
        unit = f" ({self.length_unit})" if self.length_unit else ""
        if reference == "none":
            return f"{component} position{unit}"
        name = "|displacement|" if component == "mag" else f"{component} displacement"
        return f"{name}{unit}"

    # ---- plots ---------------------------------------------------------------
    def plot_snapshot(self, at, unit="step", component="z", reference="initial",
                      x_axis="index", ax=None, save_path=None, show=True, **plot_kw):
        """Displacement of every grain at one (or several) instants.

        at     : step / time, or a list of them to overlay.
        x_axis : 'index'    -> grain number (0 = driven grain)
                 'position' -> grain x coordinate at that instant (chain shape;
                               with component='z', reference='none' this is the
                               literal x-z configuration of the chain).
        """
        ats = [at] if np.isscalar(at) else list(at)
        steps = [self.to_step(a, unit) for a in ats]
        disp = self.displacement(component, reference)
        fig, ax = _get_ax(ax, figsize=(8, 3))
        colors = _colors(len(steps))
        for c, s in zip(colors, steps):
            xs = np.arange(self.n_grain) if x_axis == "index" else self.coordinate("x")[s]
            label = f"t = {self.times[s]:.4g} s" if self.dt is not None else f"step {s}"
            ax.plot(xs, disp[s], marker="o", ms=3, lw=1.5, color=c, label=label, **plot_kw)
        ax.set_xlabel("Grain index (0 = driven grain)" if x_axis == "index"
                      else f"x position ({self.length_unit})")
        if x_axis == "index":
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylabel(self._ylabel(component, reference))
        ax.set_title(f"{self.name} - grain {component} displacement snapshot".strip(" -"))
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=8)
        return _finish(fig, ax, save_path, show)

    def plot_time_history(self, grains=None, start=None, end=None, unit="step",
                          component="z", reference="initial", ax=None,
                          save_path=None, show=True, **plot_kw):
        """Displacement vs time of selected grains over the window [start, end].

        grains : list of drive-relative indices (0 = driven); None = all grains.
        """
        sl = self.window(start, end, unit)
        gs = self._grains(grains)
        disp = self.displacement(component, reference)
        t = self.times[sl]
        fig, ax = _get_ax(ax, figsize=(8, 3))
        for c, g in zip(_colors(len(gs)), gs):
            ax.plot(t, disp[sl, g], color=c, lw=1.5, label=f"grain {g}", **plot_kw)
        ax.set_xlabel(self.time_label)
        ax.set_ylabel(self._ylabel(component, reference))
        ax.set_title(f"{self.name} - {component} displacement history".strip(" -"))
        ax.grid(True, alpha=0.4)
        if len(gs) <= 12:
            ax.legend(fontsize=8, ncol=min(len(gs), 4))
        return _finish(fig, ax, save_path, show)

    def plot_spacetime(self, start=None, end=None, unit="step", component="z",
                       reference="initial", cmap="RdBu_r", ax=None,
                       save_path=None, show=True):
        """Space-time map (grain index vs time) of the displacement over a window."""
        sl = self.window(start, end, unit)
        disp = self.displacement(component, reference)[sl]
        t = self.times[sl]
        fig, ax = _get_ax(ax, figsize=(8, 4))
        lim = np.nanmax(np.abs(disp)) or 1.0
        im = ax.pcolormesh(t, np.arange(self.n_grain), disp.T, cmap=cmap,
                           vmin=-lim, vmax=lim, shading="nearest")
        fig.colorbar(im, ax=ax, label=self._ylabel(component, reference))
        ax.set_xlabel(self.time_label)
        ax.set_ylabel("Grain index (0 = driven grain)")
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(f"{self.name} - {component} displacement space-time".strip(" -"))
        return _finish(fig, ax, save_path, show)


# ----------------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------------
def _get_ax(ax, figsize):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    return fig, ax


def _colors(n):
    cmap = plt.get_cmap("viridis")
    return [cmap(v) for v in np.linspace(0, 0.9, max(n, 1))]


def _finish(fig, ax, save_path, show):
    fig.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        fig.savefig(save_path, dpi=150)
    if show:
        plt.show()
    return fig, ax


# ----------------------------------------------------------------------------
# Command line / script entry point
# ----------------------------------------------------------------------------
# Defaults -- edit here or override on the command line.
DEFAULTS = dict(
    data_dir=os.path.dirname(os.path.abspath(__file__)),  # replace with database folder
    file="positions_406.dat",
    n_grain="auto",          # 'auto' = detect from data, or an int, e.g. 19
    drive_end="first",       # 'first' = grain 0 in file is driven, 'last' = end grain
    frame_rate=240.0,        # experiment: frames per second; set None and use dt_out for sims
    dt_out=None,             # simulation: dt * time_scale * nTout (s)
    length_unit="mm",
    component="z",
    reference="initial",
)


def _build_parser():
    p = argparse.ArgumentParser(description="Plot grain-chain displacement from a positions .dat file")
    p.add_argument("--data-dir", default=DEFAULTS["data_dir"])
    p.add_argument("--file", default=DEFAULTS["file"])
    p.add_argument("--n-grain", default=DEFAULTS["n_grain"],
                   help="number of grains per step, or 'auto' to detect (default)")
    p.add_argument("--drive-end", default=DEFAULTS["drive_end"], help="'first' (grain 0) or 'last'")
    p.add_argument("--frame-rate", type=float, default=DEFAULTS["frame_rate"])
    p.add_argument("--dt-out", type=float, default=DEFAULTS["dt_out"])
    p.add_argument("--length-unit", default=DEFAULTS["length_unit"])
    p.add_argument("--component", default=DEFAULTS["component"], choices=["x", "y", "z", "mag"])
    p.add_argument("--reference", default=DEFAULTS["reference"],
                   help="'initial', 'mean', 'none' or a step index")
    p.add_argument("--unit", default="step", choices=["step", "s"],
                   help="unit of --snapshot / --window values")
    p.add_argument("--snapshot", type=float, nargs="+", default=None,
                   help="instant(s) for the snapshot plot")
    p.add_argument("--window", type=float, nargs=2, default=None, metavar=("START", "END"),
                   help="time window for the history plot (default: whole record)")
    p.add_argument("--grains", type=int, nargs="+", default=None,
                   help="drive-relative grain indices for the history plot (default: all)")
    p.add_argument("--spacetime", action="store_true", help="also plot the space-time map")
    p.add_argument("--list", action="store_true", help="list positions_*.dat in data-dir and exit")
    p.add_argument("--save-dir", default=None, help="save figures here (PNG)")
    p.add_argument("--no-show", action="store_true", help="do not open figure windows")
    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)
    if args.list:
        for f in find_data_files(args.data_dir):
            print(f)
        return

    dt_kw = {"dt_out": args.dt_out} if args.dt_out is not None else {"frame_rate": args.frame_rate}
    chain = ChainData.from_file(args.data_dir, args.file, args.n_grain,
                                drive_end=args.drive_end, length_unit=args.length_unit, **dt_kw)
    print(chain.summary())

    ref = int(args.reference) if args.reference.lstrip("-").isdigit() else args.reference
    show = not args.no_show
    conv = (lambda v: int(v)) if args.unit == "step" else float

    def out(tag):
        return os.path.join(args.save_dir, f"{chain.name}_{tag}.png") if args.save_dir else None

    snap = [conv(v) for v in args.snapshot] if args.snapshot else [chain.n_steps // 2]
    chain.plot_snapshot(snap, unit=args.unit, component=args.component, reference=ref,
                        save_path=out("snapshot"), show=show)

    start, end = (conv(args.window[0]), conv(args.window[1])) if args.window else (None, None)
    chain.plot_time_history(args.grains, start, end, unit=args.unit, component=args.component,
                            reference=ref, save_path=out("history"), show=show)
    if args.spacetime:
        chain.plot_spacetime(start, end, unit=args.unit, component=args.component,
                             reference=ref, save_path=out("spacetime"), show=show)


if __name__ == "__main__":
    main()
