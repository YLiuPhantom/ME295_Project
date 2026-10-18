# Running Log — ME295 Project

**Breakdown of the continuum limit in loosely coupled ring chains** (see `README.md` for the research questions and approach).

This log tracks progress at the project level. Each top-level item is a project task; coding details are kept as subsections under the task they serve.

_Last updated: 2026-10-18_

---

## Status overview

| # | Project task | Status |
|---|---|---|
| 0 | Proposal, literature review, repository setup | In progress |
| 1 | Environment: LSDEM build, regression test, determinism test | Not started |
| 2 | Post-processing toolkit (data loading, basic plots, analysis) | In progress |
| 3 | Phase 1 — Map the transition: two-route λ/a sweep and collapse test (RQ1) | Not started |
| 4 | Phase 2 — Laminar detection and intermittency classification (RQ2) | Not started |
| 5 | Phase 3 — Chaos: twin-trajectory and comoving Lyapunov exponents (RQ3) | Not started |
| 6 | Phase 4 — Reduced point-mass model and bifurcation continuation (RQ4) | Not started |
| 7 | Write-up | Not started |

---

## 0. Proposal, literature review, repository setup

- [x] Draft project README: overview, research questions RQ1–RQ4, approach, roadmap, references.
- [ ] Finish the pre-proposal and the related-work note (`docs/`).
- [ ] Set up the repository layout described in the README (`lsdem/`, `analysis/`, `model/`, `data/`, `docs/`, `notebooks/`).
- [ ] Add `requirements.txt` and fill in the README TODOs (LSDEM build, baseline commands, license, contact).

## 1. Environment: LSDEM build and checks

- [ ] Build LSDEM; record compiler and MPI/OpenMP versions.
- [ ] Reproduce one baseline case from the original manuscript and compare with the reference output.
- [ ] Zero-perturbation twin test: two identical runs must stay bitwise identical (required before any Lyapunov work).
- [ ] Make each run record commit hash, input deck, and seed.

## 2. Post-processing toolkit

Goal: one module that loads any chain output (simulation or experiment) and produces the plots and order parameters needed in Phases 1–3.

### 2.1 Data loading — done (first version)
- [x] `chain_dynamics.py`: parse `positions_*.dat` (rows = grains, `n_grain` rows per output step, columns x y z) into an array `(n_steps, n_grain, 3)`.
- [x] Settable inputs: `data_dir`, `filename`, `n_grain`, `drive_end` (`"first"` / `"last"`), time base (`frame_rate` or `dt_out`), `length_unit`, `length_scale`.
- [x] Automatic `n_grain` detection (`"auto"`), with manual override and a warning on mismatch.
  - First version (global minimum over lag) picked multiples of n for perfectly periodic motion (n = 3, 5). Replaced with "smallest sharp dip vs. neighbouring lags, confirmed at 2k".
  - Passes synthetic cases (n = 2…5000, straight and curved chains, truncated last step) and the real file (dip score 0.35 vs ≈20 at neighbouring lags).
- [x] Grain 0 is always the driven grain in analysis functions; `raw_index()` maps back to file order.
- [x] Check for repeated frames.
- [ ] Load rotations / forces / intergrain_forces files with the same loader (4- and 5-column versions).
- [ ] Cache parsed arrays (`.npy` / `.npz`) for faster reloads of large files.

### 2.2 Basic plots — done (first version)
- [x] Displacement of every grain at chosen time(s) (x-axis = grain index or x position).
- [x] Displacement of chosen grains over a time window (steps or seconds).
- [x] Space-time map and chain-shape figure; verified visually on `positions_406.dat`.
- [x] Command-line interface (`--help`).
- [ ] Amplitude-vs-position (envelope) plot and comparison with theory.
- [ ] Animation / chain-shape movie over a time window.

### 2.3 Analysis functions — to do
- [ ] Port analysis from the reference script into the new module: velocity / KE, spectra (PSD, MFER, spectral entropy), wavelets, intermittency measures.
- [ ] Add order parameters needed for Phase 1: slack-link fraction φ<sub>s</sub>, snap-event rate.
- [ ] Save processed data to CSV.

### 2.4 Data management and testing
- [ ] Define naming / folder layout for the database (e.g. `<case>/Output/positions_*.dat`) and add a batch loop over cases.
- [ ] Unit tests (pytest): loader, `n_grain` detection, drive-end flip, window handling.
- [x] Usage notes: `README_chain_dynamics.md` (examples run against `positions_406.dat`).

### 2.5 Open questions / bottlenecks
- **Units and time base** not yet confirmed per source (experiment: mm @ 240 Hz assumed; simulation: voxel, `dt*time_scale*nTout`).
- **Repeated frames** in experiment files (42 in `positions_406.dat`): keep, mask, or interpolate?
- **Drive end per case** should come from case metadata; currently set by hand. In the example, grain 0 has the largest z amplitude, consistent with `drive_end="first"`.
- **Load time:** `pandas.read_csv` reads the full text file each call; may be slow for full simulation outputs → caching (2.1).
- **`n_grain` auto-detection limits:** assumes grains move less than about half a spacing per output step; fails (with a clear error) for single-step or very coarse output — pass `n_grain` manually.

## 3. Phase 1 — Map the transition (RQ1)

- [ ] Route A: sweep drive parameters at fixed geometry.
- [ ] Route B: sweep ring geometry (spacing a) at fixed drive.
- [ ] Compute order parameters (MFER, φ<sub>s</sub>, snap rate, spectral content) vs λ/a; test whether both routes collapse.
- [ ] Fine sweep near integer λ/a to test commensurability locking.

## 4. Phase 2 — Intermittency (RQ2)

- [ ] Laminar-window detection: stroboscopic return distance + wavelet-energy criterion.
- [ ] Measure P(τ) and ⟨τ⟩(ε); compare with intermittency types I, II, III, V, on-off.
- [ ] Apply the same pipeline to experimental records.

## 5. Phase 3 — Chaos (RQ3)

- [ ] Maximal Lyapunov exponent (twin trajectory + Benettin renormalization) — requires the determinism test (task 1).
- [ ] Finite-time exponents in laminar vs burst windows; lifetime statistics to rule out transient chaos.
- [ ] Comoving Lyapunov exponents Λ(v): absolute vs convective chaos.
- [ ] Experiments: 0-1 test and recurrence quantification.

## 6. Phase 4 — Reduced model (RQ4)

- [ ] Point-mass chain with unilateral links and impact law; event-driven integrator.
- [ ] Saltation matrices for exact variational dynamics.
- [ ] Continue periodic waves in λ/a; locate the bifurcation at λ<sub>c</sub>/a.

## 7. Write-up

- [ ] Progress reports, figures, final report.

---

## Dated log

### 2026-10-18
- Inspected `positions_406.dat`: 22 800 rows × 3 cols → 19 grains × 1200 steps; y ≡ 0 (x-z plane motion); x increases with grain index in every frame; mean spacing ≈ 17 units.
- Reviewed the reference post-processing script for conventions (displacement = coordinate minus mean/initial; grains labelled `n_grain-1-g` when driven from the end).
- Wrote `chain_dynamics.py` (`ChainData` class + CLI) and `README_chain_dynamics.md` (task 2.1–2.2).
- Improved `n_grain` auto-detection after failures on perfectly periodic synthetic chains.
- Expanded the project `README.md` (overview, RQs, approach, roadmap).
- Added `.gitignore` (`*.dat`, `*.docx`, `YL_postprocess.py`, `__pycache__/`); pushed running log to GitHub.
