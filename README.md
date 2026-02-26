# Time-varying directed interaction estimation

MATLAB reference implementation of the **Sliding-Window Prediction Correlation (SWpC)** framework for estimating the **time-varying strength** (prediction-correlation) and **duration** (how long directed influence persists across windows) of directed interactions between ROI time series (e.g., fMRI BOLD).

If you use this code, please cite the SWpC paper:  
**Paper/Preprint:** Xu N, Zhang X, Pan WJ, Smith JL, Schumacher EH, Allen JW, Calhoun VD, Keilholz SD. *Time-Varying Directed Interactions in Functional Brain Networks: Modeling and Validation.* bioRxiv (2026). https://doi.org/10.64898/2026.02.15.706049  
**Code:** https://github.com/inspirelab-site/swpc

---

## Contents

- `main_script2run_swpc.m` — demo driver script (load data → run SWpC → save outputs)
- `data.mat` — demo input data
- `swpcv0323/`
  - `sliding_pcorr_window.m` — sliding-window wrapper
  - `pcorr_x2y.m` — core estimator (fit + order selection + correlation)

---

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox (`aicbic`)
- Optimization Toolbox only if using positivity constraint (`lsqlin` via `Constraint='p'`)

---

## Quick start

In MATLAB (repo root):

```matlab
run('main_script2run_swpc.m')
```
---

## Using your own data

SWpC is called on **two vectors per direction**:

- `tc1` (source): `(T x 1)`
- `tc2` (target): `(T x 1)`

To use more ROIs, construct input as `(nroi x T)` and loop over ROI pairs (as in the demo script).

---


## Key parameters (in `main_script2run_swpc.m`)

- `TR` (sec)  
  Sampling interval.

- `windowsize_duration` (sec)  
  Window duration in seconds. Window length in samples is computed from `TR` and forced even:
  - `windowsize = floor(windowsize_duration/TR/2)*2`

- `overlap` (samples)  
  Window overlap. The demo uses near full overlap:
  - `overlap = windowsize - 1`

- `maxDur` (sec) and `maxLag` (samples)  
  Maximum lag/order search range:
  - `maxLag = floor(maxDur/TR)`

- `Constraint`  
  Regression mode used to estimate the impulse response:
  - `'p'`  : positivity-constrained least squares (uses `lsqlin`)
  - `'all'`: unconstrained least squares (uses `lsqr`)

---
## Method overview

Within each sliding window, SWpC estimates directed influence **tc1 → tc2** by:

1. Constructing a **lag-embedded predictor matrix** from the source window (`tc1`) up to candidate model orders (`Nh = 1..maxLag`).
2. Fitting a **finite impulse response (FIR) filter** `h` (lag weights) to predict the target window (`tc2`) from the lagged source.
3. Selecting the model order `Nh` per window using information criteria (**BIC** and **AIC/AICc**).
4. Quantifying directed interaction strength in that window as the **prediction correlation** between the observed target (`tc2`) and the model-predicted target (`tc2_hat`).

Running this window-by-window yields a **time series of directed strength** for each ROI pair and direction.

---

## Strength vs. duration

How to interpret SWpC outputs:

- **Strength**: the SWpC value per window (e.g., `corr_bic` or `corr_aic`) is the instantaneous directed interaction strength at that window.
- **Duration**: duration is not returned as a single primitive value; it is derived from the **temporal persistence** of directed strength across windows. For example, you can compute duration as:
  - the total number of windows where `corr_*` exceeds a threshold, and/or
  - the length of contiguous above-threshold segments (then convert windows → seconds using `TR` and the window step).

---

## AIC vs AICc

This implementation computes both **AIC** and **AICc** (small-sample corrected AIC).

The returned AIC-based choice (`Nh_aic`, `corr_aic`) should be interpreted as:

- **AICc-selected** in short-window regimes (when the sample size per window is not large relative to the number of fitted parameters),
- otherwise **AIC-selected**.

AICc is recommended for short windows because it penalizes model complexity more strongly and tends to yield more stable (less overfit) order selection.

---

## Outputs

### What `sliding_pcorr_window` returns (per direction)

Calling:

```matlab
[corr_bic, Nh_bic, corr_aic, Nh_aic, corr_std, corr0, H, nmaps] = ...
    sliding_pcorr_window(tc1, tc2, windowsize, overlap, maxLag, Constraint);
```
Returns (typically) `1 × nWindows` vectors (unless otherwise noted):

- `corr_bic` — prediction correlation using **BIC-selected** order (per window)
- `Nh_bic` — **BIC-selected** order `Nh` (per window)
- `corr_aic` — prediction correlation using **AIC/AICc-selected** order (per window)
- `Nh_aic` — **AIC/AICc-selected** order `Nh` (per window)
- `corr_std` — baseline correlation (see `pcorr_x2y.m` for exact definition used)
- `corr0` — prediction correlation using order `Nh = 1`
- `H` — estimated impulse responses per window (stacked across windows)
- `nmaps` — number of windows

### What `main_script2run_swpc.m` saves (per subject)

The demo script saves one file per subject:

- `./output/subj<k>.mat`

Key saved arrays are shaped:

- `(nroi x nroi x nWindows)` with index order: **`[sourceROI, targetROI, window]`**

Saved variables include:

- `slpcorr_BIC`, `slNh_BIC`
- `slpcorr_AIC`, `slNh_AIC` *(AIC/AICc per the rule above)*
- `slpcorr_STD`
- `slpcorr1` *(stores `corr0`)*

Both directions are computed and stored (**i→j** and **j→i**).
