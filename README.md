# EigenDataAnalysis.jl

**EigenDataAnalysis.jl** is a lightweight, high-performance Julia package for Operational Modal Analysis (OMA) and Stochastic Subspace Identification (SSI). It provides tools to extract modal parameters (natural frequencies, damping ratios, and mode shapes) and system state-space models directly from output-only response data.

---

## Key Features

* **3 Core SSI Algorithms:**
  * `SSICOV`: Covariance-driven SSI using Block-Hankel covariance matrices.
  * `PCSSI`: Principal Component SSI using Block-Toeplitz data projections.
  * `SSIDATA`: Data-driven SSI using orthogonal projection via LQ decomposition.
* **Stabilization Diagrams:** Automated pole convergence tracking across varying model orders (`:fully_stable`, `:freq_damp`, `:freq_mac`, `:freq`, `:unstable`).
* **Modal Validation:** Built-in Modal Assurance Criterion (`mac`) and modal damping ratio estimation.
* **Response Reconstruction:** Time-domain signal reconstruction and quantitative error estimation (`reconstruction_error`).
* **Unstable Pole Filtering:** Built-in detection and optional filtering (`unstable_filt=true`) for non-physical unstable poles ($Re(\omega) > 0$).

---

## Installation

Until published to the General Registry, load the package locally in Julia:

```julia
using Pkg
Pkg.activate(".")
using EigenDataAnalysis

```

---

## Quick Start

```julia
using EigenDataAnalysis
using LinearAlgebra

# 1. Prepare output time-series data X (ny × N) and sampling interval dt
ny, N = 4, 2000     # 4 sensor channels, 2000 time points
dt = 0.01           # 100 Hz sampling rate
X = randn(ny, N)    # Measured response matrix

# 2. Fit a system identification model using SSI-COV
model = fit(SSICOV, X, dt; max_lag=20, r=6, unstable_filt=true)

# 3. Extract modal properties
ω = eigenvalues(model)          # Continuous-time complex poles
f_hz = imag.(ω) ./ (2π)         # Natural frequencies (Hz)
zetas = damping_ratios(model)   # Damping ratios
Φ = modes(model)                # Mode shape vectors (ny × r)
S = singular_values(model)      # Singular value spectrum

# 4. Reconstruct response signal and evaluate error
X_rec, t_idx = reconstruct(model)
err = reconstruction_error(model, X)
println("Relative Reconstruction Error: ", round(err * 100, digits=2), "%")

# 5. Compute stabilization diagram data across model orders
points, recon_errors, orders = compute_stabilization_diagram(
    SSICOV, X, dt;
    orders=2:2:20,
    max_lag=20,
    tol_f=0.01,
    tol_zeta=0.05,
    tol_mac=0.95
)

```

---

## API Reference

### Methods

| Method Type | Description |
| --- | --- |
| `SSICOV` | Covariance-Driven SSI (Block-Hankel matrix built from response covariances) |
| `PCSSI` | Principal Component SSI (Block-Toeplitz projection of past and future data) |
| `SSIDATA` | Data-Driven SSI (LQ decomposition of block-Hankel response matrix) |

### Key Functions

* **`fit(Method, X, dt; max_lag=20, r=10, unstable_filt=false)`**
Estimates the state-space model and modal parameters from response matrix $X$ ($n_y \times N$) sampled at time interval $dt$.
* **`compute_stabilization_diagram(Method, X, dt; orders, tol_f=0.01, tol_zeta=0.05, tol_mac=0.95, kwargs...)`**
Computes pole stability tracking across specified system orders (`orders`). Returns stabilization points, reconstruction errors, and evaluated orders.
* **`reconstruct(model::SSIResult)`**
Returns `(output_sequence, time_indices)` representing the reconstructed response signal and valid time span.
* **`reconstruction_error(model::SSIResult, X_original::AbstractMatrix)`**
Computes the relative normalized error $\Vert{} X_{\text{orig}} - X_{\text{rec}} \Vert{} / \Vert{} X_{\text{orig}} \Vert{}$ aligned with the model output timeframe.
* **`mac(ϕ1, ϕ2)`**
Computes the Modal Assurance Criterion scalar ($0 \le \text{MAC} \le 1$) between two mode vectors.
* **`damping_ratios(model::SSIResult)`**
Calculates the viscous damping ratios ($\zeta = -Re(\omega)/\vert{}\omega\vert{}$) for all poles.

### Accessor Functions

* **`singular_values(res)`**: Returns singular value vector $S$.
* **`system_matrices(res)`**: Returns a named tuple `(A = res.A, C = res.C)`.
* **`modes(res)`**: Returns mode shape matrix $\Phi$.
* **`eigenvalues(res)`**: Returns continuous-time complex eigenvalues $\omega$.

---

## Repository Structure

```text
EigenDataAnalysis/
├── src/
│   └── EigenDataAnalysis.jl
├── notebooks/
│   └── ssi_testing.jl
├── README.md
├── LICENSE
└── .gitattributes

```

```

```
