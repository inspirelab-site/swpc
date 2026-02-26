# SWpC (MATLAB): Sliding-Window Prediction Correlation

MATLAB reference implementation of **Sliding-Window Prediction Correlation (SWpC)** for estimating **time-varying directed interactions** between ROI time series (e.g., fMRI BOLD).

---

## Using your own data

SWpC is called on **two vectors per direction**:

- `tc1` (source): `(T x 1)`
- `tc2` (target): `(T x 1)`

To use more ROIs, construct `bold` as `(nroi x T)` and loop over ROI pairs (as in the demo script).

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

Within each sliding window, SWpC:

1. **Lag-embeds** the source window to form a design matrix of delayed versions of the source time series.
2. **Fits an impulse response / filter** `h` such that the lag-embedded source predicts the target window.
3. **Selects model order** (`Nh`, number of lags) per window using information criteria.
4. **Computes prediction correlation** between the observed target window and the predicted target window.

Repeating across windows yields a time-resolved directed interaction estimate.

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
