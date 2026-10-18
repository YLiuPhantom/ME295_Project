# chain_dynamics.py

Reads grain-chain position files (`positions_*.dat`) from LSDEM simulations or tracked experiments and plots grain displacements:

- the displacement of **every grain at a chosen time** (snapshot), and
- the displacement of **chosen grains over a chosen time window** (time history),

plus a space-time map of the whole chain.

Requirements: Python ≥ 3.9, `numpy`, `pandas`, `matplotlib`.

---

## 1. Input data format

Each `.dat` file is plain text with whitespace-separated columns `x y z`. Every row is one grain. Each output time step is a block of `n_grain` consecutive rows:

```
x_0 y_0 z_0      <- step 0, grain 0
x_1 y_1 z_1      <- step 0, grain 1
...
x_N-1 ...        <- step 0, grain n_grain-1
x_0 y_0 z_0      <- step 1, grain 0
...
```

The loader reshapes this into an array of shape `(n_steps, n_grain, 3)`. If the file ends partway through a step (e.g. a simulation that was stopped), the incomplete last step is dropped with a warning.

Example: `positions_406.dat` has 22 800 rows, i.e. 19 grains × 1200 steps.

---

## 2. Key conventions

| Item | Convention |
|---|---|
| **Grain numbering** | In all analysis and plot functions, **grain 0 = the driven grain** and grain `n_grain-1` = the far end. If the drive is at the last grain in the file, the grain order is flipped internally. `chain.raw_index(g)` returns the row position in the file. |
| **Displacement reference** | `reference="initial"` (default): position minus position at step 0. `"mean"`: minus each grain's time mean (as in `YL_postprocess.py`). An integer `k`: minus position at step `k`. `"none"`: raw coordinates. |
| **Components** | `"x"`, `"y"`, `"z"`, or `"mag"` (length of the displacement vector). |
| **Time units** | Times can be given in output steps (`unit="step"`, default) or seconds (`unit="s"`, needs a time base). |
| **Length units** | Data is not rescaled unless you set `length_scale`. `length_unit` only sets the axis label. |

---

## 3. Variables you can change

| Variable | Meaning | Default |
|---|---|---|
| `data_dir` | Folder that holds the `.dat` files (your database, e.g. `D:\LSDEM_DataCenter\test\<case>\Output`) | folder of the script |
| `filename` / `--file` | File inside `data_dir` | `positions_406.dat` |
| `n_grain` | Number of grains in the chain. `"auto"` detects it from the data, or give an int (e.g. `19`). If your value disagrees with the detected one, you get a warning and your value is used. | `"auto"` |
| `drive_end` | Which end is driven: `"first"` (grain 0 in the file) or `"last"` (end grain). Also accepts `0`/`-1`, `"start"`/`"end"`. | `"first"` |
| `frame_rate` | Samples per second (experiments, e.g. 240 Hz) | 240 (CLI); `None` (Python) |
| `dt_out` | Seconds between output steps (simulations: `dt * time_scale * nTout`). Use **either** `frame_rate` **or** `dt_out`. If neither is set, time is in output steps. | `None` |
| `length_unit` | Axis label for lengths (`"mm"`, `"voxel"`, `"m"`, …) | `"mm"` |
| `length_scale` | Multiply raw coordinates by this (e.g. `4e-4` for voxel → m) | `1.0` |

There are three ways to set them:

1. **On the command line** (see §5). Nothing in the file needs editing.
2. **Edit the `DEFAULTS` dictionary** near the bottom of `chain_dynamics.py`. This sets the defaults used when you run the script with no arguments.
3. **In your own Python script**, as arguments to `ChainData.from_file(...)` (see §6).

---

## 4. Functions and methods

### Module-level functions

| Function | What it does |
|---|---|
| `load_positions(file_path, n_grain="auto")` | Reads one `.dat` file and returns an array of shape `(n_steps, n_grain, 3)` in file order. |
| `detect_n_grain(rows)` | Detects the number of grains per step from the raw rows (see §7). |
| `find_data_files(data_dir, pattern="positions_*.dat")` | Returns a sorted list of matching files in a folder. |
| `normalize_drive_end(drive_end)` | Converts the different ways of writing the drive end to `"first"` or `"last"`. |

### `ChainData` class

Create it with `ChainData.from_file(data_dir, filename, n_grain="auto", drive_end="first", frame_rate=..., ...)` or `ChainData(positions_array, ...)`.

**Attributes**

| Name | Description |
|---|---|
| `positions` | `(n_steps, n_grain, 3)`, ordered so grain 0 is the driven grain |
| `raw_positions` | Same data in file order |
| `n_steps`, `n_grain` | Array sizes |
| `times` | Time of each step (s), or step index if no time base was given |

**Data methods**

| Method | Returns |
|---|---|
| `displacement(component="z", reference="initial")` | `(n_steps, n_grain)` displacement array |
| `coordinate(component="z")` | `(n_steps, n_grain)` absolute coordinate |
| `to_step(value, unit="step"/"s")` | Converts a time or step to a valid step index (negative values count from the end) |
| `window(start, end, unit)` | Step slice `[start, end]`, inclusive |
| `raw_index(g)` | Converts a grain index (0 = driven grain) to its row index in the file |
| `repeated_frames()` | Steps identical to the step before (e.g. dropped frames in tracking) |
| `summary()` | One-line description: size, drive end, duration, repeated frames |

**Plot methods.** All plot methods accept `ax=` (draw on an existing axes), `save_path=` (save a PNG), and `show=` (open a window or not). They return `(fig, ax)`.

| Method | Plot |
|---|---|
| `plot_snapshot(at, unit="step", component="z", reference="initial", x_axis="index")` | Displacement of all grains at time `at`. Pass a list to overlay several times. Use `x_axis="position"` to plot against each grain's x coordinate. With `component="z", reference="none", x_axis="position"` you get the actual chain shape in the x-z plane. |
| `plot_time_history(grains=None, start=None, end=None, unit="step", component="z", reference="initial")` | Displacement vs time for the chosen grains (`None` = all) between `start` and `end`. |
| `plot_spacetime(start=None, end=None, unit="step", component="z", reference="initial")` | Colour map of displacement, grain index vs time. |

---

## 5. Running from the command line

```powershell
# Defaults (example file in this folder, n_grain auto, drive at grain 0, 240 Hz)
python chain_dynamics.py

# Real database, z snapshots at 1.25 s and 2.5 s, history of grains 0/6/12/18 for 0–2.5 s,
# plus the space-time map, figures saved to .\figs without opening windows
python chain_dynamics.py --data-dir "D:\LSDEM_DataCenter\test\exp_disconnected_stable\Output" `
    --file positions_406.dat --drive-end first --frame-rate 240 `
    --unit s --snapshot 1.25 2.5 --window 0 2.5 --grains 0 6 12 18 `
    --spacetime --save-dir figs --no-show

# Simulation output: manual grain count, drive at the end grain, time base from dt*nTout
python chain_dynamics.py --data-dir "D:\LSDEM_DataCenter\test\case\Output" --file positions_406-full.dat `
    --n-grain 21 --drive-end last --dt-out 4.36916e-4 --length-unit voxel

# List the positions files available in a folder
python chain_dynamics.py --data-dir "D:\LSDEM_DataCenter\test\case\Output" --list
```

| Option | Meaning |
|---|---|
| `--data-dir`, `--file` | Data folder and file name |
| `--n-grain` | Integer, or `auto` (default) |
| `--drive-end` | `first` or `last` |
| `--frame-rate` / `--dt-out` | Time base. `--dt-out` takes priority when given. |
| `--length-unit` | Axis label for lengths |
| `--component` | `x`, `y`, `z`, or `mag` |
| `--reference` | `initial`, `mean`, `none`, or a step index |
| `--unit` | `step` or `s`, used for `--snapshot` and `--window` |
| `--snapshot T1 [T2 ...]` | Time(s) for the snapshot plot (default: the middle step) |
| `--window START END` | Window for the history and space-time plots (default: the whole record) |
| `--grains G1 [G2 ...]` | Grains for the history plot, 0 = driven grain (default: all) |
| `--spacetime` | Also draw the space-time map |
| `--save-dir` | Save PNGs as `<file>_snapshot.png`, `<file>_history.png`, `<file>_spacetime.png` |
| `--no-show` | Don't open figure windows (for batch runs) |
| `--list` | List `positions_*.dat` files in `--data-dir` and exit |

Run `python chain_dynamics.py --help` for the full list.

---

## 6. Using it from Python

```python
from chain_dynamics import ChainData, find_data_files

data_dir = r"D:\LSDEM_DataCenter\test\exp_disconnected_stable\Output"

chain = ChainData.from_file(data_dir, "positions_406.dat",
                            n_grain="auto",      # or 19
                            drive_end="first",   # or "last"
                            frame_rate=240.0)    # or dt_out=dt*time_scale*nTout
print(chain.summary())

# Displacement of every grain at t = 1.25 s and 2.5 s
chain.plot_snapshot(at=[1.25, 2.5], unit="s")

# Chain shape in the x-z plane at step 300
chain.plot_snapshot(at=300, reference="none", x_axis="position")

# z displacement of grains 0, 6, 12 between steps 0 and 600
chain.plot_time_history(grains=[0, 6, 12], start=0, end=600)

# Raw arrays for your own analysis
z = chain.displacement("z", reference="mean")   # shape (n_steps, n_grain)
t = chain.times

# Loop over every case file in a folder
import os
for path in find_data_files(data_dir):
    c = ChainData.from_file(data_dir, os.path.basename(path), frame_rate=240)
    c.plot_time_history(grains=[0], show=False, save_path=f"figs/{c.name}_g0.png")
```

---

## 7. How `n_grain` auto-detection works

Rows `i` and `i + n_grain` are the **same grain one step later**, so they are close together (distance = motion per step). Rows at any other lag are **different grains** (distance ≈ grain spacing). The code computes the median row-to-row distance at each lag. It picks the smallest lag `k ≥ 2` whose distance is a sharp drop compared with lags `k-1` and `k+1`, and that repeats at `2k`.

- On `positions_406.dat` the distance at lag 19 is 0.35, against about 20 at lags 18 and 20.
- Tested on synthetic chains of 2–5000 grains; detection takes under 1 s.
- **Limitation:** it needs grains to move less than about half a grain spacing between saved steps, and at least two steps in the file. Otherwise it stops with a clear error and you pass `n_grain` by hand.

---

## 8. Notes and known limitations

- Units and time base are **not** read from the file. Set `frame_rate`/`dt_out` and `length_unit` for each data source.
- `drive_end` is set by hand for each file. In the example, grain 0 has the largest z amplitude, which fits `drive_end="first"`.
- Experiment files can contain repeated frames (42 in the example). They are counted by `repeated_frames()` but not removed.
- The whole text file is parsed on every load, which may be slow for very large simulation outputs. Caching is on the to-do list in `RUNNING_LOG.md`.
- Planned additions (rotations/forces loaders, spectra, wavelets, batch processing, tests) are tracked in `RUNNING_LOG.md`.
