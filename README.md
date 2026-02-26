# SWpC (MATLAB): Sliding-Window Prediction Correlation

This repository provides a MATLAB reference implementation of **Sliding-Window Prediction Correlation (SWpC)** for estimating **time-varying directed interactions** between ROI time series.

SWpC estimates *directional influence* by predicting a target time course from a **lag-embedded** source time course within each sliding window, then computing the **prediction correlation**. The model order (lag length) can be selected per window using **AIC** or **BIC**.

Repo link: https://github.com/inspirelab-site/swpc

---

## What’s inside

- `main_script2run_swpc.m`  
  End-to-end example: loads the demo data, runs SWpC for ROI pairs, and saves per-subject outputs.

- `data.mat`  
  Example input data used by the main script.

- `swpcv0323/`  
  Core functions:
  - `sliding_pcorr_window.m` — runs SWpC across sliding windows
  - `pcorr_x2y.m` — fits the lagged prediction model, computes prediction correlation, performs AIC/BIC order selection

- `output/` *(created when you run the script)*  
  Per-subject results saved as `.mat` files (e.g., `subj1.mat`).

---

## Requirements

- MATLAB (recent versions recommended)
- Toolboxes (depending on settings):
  - **Statistics and Machine Learning Toolbox** (for `aicbic`)
  - **Optimization Toolbox** (for constrained least squares via `lsqlin`, used when `Constraint='p'`)

> If you do not have Optimization Toolbox, set `Constraint='all'` in the main script to avoid `lsqlin`.

---

## Quick start

1) Clone the repo and open MATLAB in the repo root.

2) Run:
```matlab
run('main_script2run_swpc.m')
