# CJSSRRR_hybrid: a refined CJ--SS--RR algorithm with hybrid restarts

A MATLAB implementation of the refined CJ--SS--RR method with hybrid restart for computing eigenpairs of Hermitian matrices within a given interval.


## Overview

CJSSRRR_hybrid computes eigenpairs of Hermitian matrices that lie within a given interval `[a, b]`. It supports natural, combination, and hybrid restarts, where the default is hybrid restart.


## Function Signature

```matlab
[eigenvalues, eigenvectors, flag, MVs, info] = CJSSRRR_hybrid(A, interval, varargin)
```


## Usage Examples

#### Basic Usage

```matlab
% Create a test matrix
% The discretization of a differential operator
% on a regular two-dimensional lattice
L = gallery('tridiag', 50);
A = kron(speye(size(L)), L) + kron(L, speye(size(L)));

% Compute eigenvalues in the interval [1, 1.5]
% There are 119 eigenvalues of A in this interval
interval = [1, 1.5];
[E, X, flag, MVs, info] = CJSSRRR_hybrid(A, interval);

fprintf('Found %d eigenvalues.\n', length(E));
fprintf('Total iterations: %d.\n', info.iteration);
```

#### Using Structure Input

```matlab
% Custom parameters using structure
options = struct('maxItr', 50, 'moment', 8, 'eigTol', 1e-12, 'disp', 2);
[E, X, flag, MVs, info] = CJSSRRR_hybrid(A, interval, options);
```

#### Using Name-Value Pairs

```matlab
% Custom parameters using name-value pairs
[E, X, flag, MVs, info] = CJSSRRR_hybrid(A, interval, ...
    'maxItr', 50, 'moment', 8, 'eigTol', 1e-12, 'disp', 2);
```


## Input Parameters

### Required Parameters

- `A`: n-by-n Hermitian matrix.
- `interval`: 2-element vector `[a, b]` specifying the eigenvalue search interval.

### Optional Parameters (Name-Value Pairs)

Optional parameters use native inputParser structure expansion: supply a scalar structure, name-value pairs, or a mixture following the parser's native rules. Full names are case-insensitive; unknown names remain ignored and partial matching is disabled.

| Parameter | Type | Default | Description |
|:----------|:-----|:--------|:------------|
| `moment` | integer | 8 | Number of moments. |
| `maxItr` | integer | 50 | Maximum number of iterations. |
| `eigTol` | double | 1e-12 | Convergence tolerance. |
| `eigenspace` | matrix | [] | Desired eigenspace for testing purposes. |
| `spectrum` | vector | [] | Spectrum interval of A `[lambda_min, lambda_max]`. |
| `nev` | integer | 0 | Number of eigenvalues in `interval` (0 = auto-estimate). |
| `V` | matrix | [] | Starting matrix to generate the search subspace. |
| `degree` | integer | 0 | Degree of the CJ polynomial series (0 = auto-estimate). |
| `degreeD` | real scalar | 1 | Positive D parameter in degree estimation (recommended: 1–10). |
| `omega` | double | 1.5 | Scaling factor for search subspace dimension (`dim ≈ nev × omega`). |
| `disp` | integer | 0 | Diagnostic information display level (any finite real integer). |
| `restart` | string or char | `"hybrid"` | Restart format: `"natural"`, `"combine"`, or `"hybrid"`. |
| `orthMethod` | string or char | `"QR"` | Search-space orthonormalization: `"TSVD"` or `"QR"`. |
| `svdTol` | double | `100*eps` | Relative TSVD threshold based on the Frobenius norm of `S`. |
| `eps_switch` | double | `eigTol^.5` | Real scalar hybrid switching threshold in `[0, 1)`. |

### Validation and Defaults

The parser rejects incorrect types or shapes, nonfinite or complex scalar options, and noninteger values for `moment`, `maxItr`, `nev`, `degree`, and `disp`. Text options must be character row vectors or string scalars; `V` and `eigenspace` must be numeric matrices. `spectrum` must be numeric and either empty or a real vector.

Recoverable values issue a parameter-specific warning and use the documented default: nonpositive `moment`, `degreeD`, or `omega`; negative `maxItr` or `degree`; `nev` outside `[0,size(A,2)]`; `eigTol`, `svdTol`, or `eps_switch` outside `[0,1)`; and unsupported restart or orthonormalization names. `degreeD` accepts positive real values without an upper bound. `0<omega<1` warns but keeps the supplied value. `disp` has no range check or fallback.

Invalid row counts or nonfinite entries in nonempty `V` and `eigenspace` reset them to `[]`, enabling random initialization or disabling eigenspace diagnostics. A nonempty `spectrum` with fewer than two elements, nonfinite entries, or equal bounds resets to `[]` for spectrum estimation. Valid spectral vectors are reduced to their minimum and maximum.

Omitted `eps_switch` silently uses `sqrt(eigTol)` after correcting `eigTol`; an explicitly out-of-range value warns and uses the same default. Explicit zero remains valid for `maxItr`, `eigTol`, `svdTol`, and `eps_switch`.

### Parameter Details

#### Core Algorithm Parameters

- `moment`: Number of moments.
- `maxItr`: Maximum number of iterations. Set to 0 to perform the preliminary estimation steps without entering the main iteration; the eigenvalue count is estimated when `nev=0`.
- `eigTol`: Finite real scalar tolerance in `[0,1)` for relative residual norms; zero remains supported.

#### Spectral Information

- `spectrum`: Spectrum interval of A. If unknown, it will be estimated by the Lanczos method. If provided, `[lambda_min, lambda_max]` should enclose the spectrum of A.
- `nev`: Number of eigenvalues in the interval. If set to 0, it will be automatically estimated by a stochastic trace estimator based on the CJ polynomial series.

#### Subspace Construction

- `V`: Starting matrix. If empty, a Gaussian random matrix will be generated.
- `omega`: Controls the dimension of the search subspace relative to `nev`.

#### Polynomial Approximation

- `degree`: Degree of the Chebyshev--Jackson polynomial series. If set to 0, the degree will be automatically estimated.
- `degreeD`: Parameter for automatic degree estimation using the formula:

$$
d = \left\lceil {\color{red}D} \cdot \frac{ \sqrt[3]{2} \pi^2 M^{1/6} }{(b-a)^{4/3}} \right\rceil - 2.
$$

Here `M` is the number of moments and `D = degreeD`. The endpoints `a` and `b` are those of the search interval after mapping the spectrum of A to `[-1, 1]`.

#### Restart and Orthonormalization

- `restart`: `"natural"` refers to the natural restart; `"combine"` refers to the combination restart; `"hybrid"` refers to the hybrid restart.
- `eps_switch`: Finite real scalar threshold with `0 <= eps_switch < 1`. The default is computed from the validated or corrected `eigTol` as `eigTol^.5`.
- `orthMethod`: `"QR"` uses an economic QR factorization of `S`; `"TSVD"` uses an economic SVD with the truncation described below. With `moment=1`, the effective method is always `"QR"`.
- `svdTol`: Finite real scalar with `0 <= svdTol < 1`. For `svdTol>0`, retain left singular vectors satisfying $\sigma_j\geq\|S\|_F\,\mathrm{svdTol}$. Setting `svdTol=0` keeps all left singular vectors returned by the economic SVD, without truncation.

`svdTol` is used only by TSVD, and `eps_switch` only by hybrid restart.

#### Diagnostic Output

- `disp`: Any finite real integer, including negative values and values above 3, is accepted unchanged without a range warning. The usual display levels are:
  - 0: No output.
  - 1: Basic algorithm information and final results.
  - 2: Iteration-by-iteration progress, including the actual RR projection dimension.
  - 3: Detailed iteration information including subspace angles when `eigenspace` is supplied.


## Output Parameters

- `eigenvalues`: Retained refined Ritz values within the interval.
- `eigenvectors`: Corresponding refined Ritz vectors.
- `flag`: Algorithm exit flag.
  - -1: Natural restart residual growth when `eigTol > eps`.
  - 0: Converged.
  - 1: Maximum iterations reached.
  - 2: No iterations performed.
- `MVs`: Matrix-vector product count, `[actual_MVs, equivalent_MVs]`.
  - `actual_MVs`: Direct matrix-vector products with A.
  - `equivalent_MVs`: Other flops counted in the CJ approximation, converted to MVs.
- `info`: Structure containing additional information.
  - `time`: Timing information for various stages.
    - `spectrum`: Time for spectrum estimation.
    - `CJcoeffs`: Time for CJ coefficients computation.
    - `nev`: Time for eigenvalue-count estimation.
    - `subspace`: Time for subspace generation, `[MVs_time, flops_time]`.
    - `orth`: Time for full search-space orthonormalization by QR or TSVD.
    - `RR`: Time for the Rayleigh--Ritz procedure.
    - `refined`: Time for the refined procedure, `[preprocess_time, compute_time]`.
    - `residual`: Time for residual computation and associated vector updates.
    - `total`: Total time for the algorithm.
  - `spectrum`: Supplied or estimated spectrum interval of A, `[lambda_min, lambda_max]`.
  - `degree`: Degree of the CJ polynomial series.
  - `degreeD`: D parameter used in degree estimation.
  - `nev`: Supplied or estimated number of eigenvalues in the interval, not the number of converged eigenvalues.
  - `V`: Starting matrix V used for the first iteration.
  - `moment`: Number of moments used in the algorithm.
  - `dimension`: Number of columns in the initial raw search matrix `S`.
  - `iteration`: Number of iterations performed.
  - `restart`: Normalized requested restart format.
  - `orthMethod`: Effective search-space orthonormalization method.
  - `svdTol`: Relative TSVD truncation threshold.
  - `eps_switch`: Hybrid switching threshold.
  - `switch`: Hybrid switch information.
    - `itr`: Iteration at which switching occurs; 0 if no switch occurs.
    - `cnd`: 0 for no switch, 1 for the residual-threshold condition, or 2 for the residual-growth condition.
  - `itrInfo`: Structure containing iteration-wise information.
    - `dimRR`: Actual Rayleigh--Ritz projection dimension, `size(Q,2)`.
    - `residual`: Maximum relative residual norm of the retained refined Ritz pairs.
    - `nev`: Number of retained refined Ritz values.
    - `condition`: Condition numbers of the retained refined Ritz vectors.
    - `insideRes`: Maximum relative residual norm of refined Ritz pairs whose values lie within the interval, before removing copies.
    - `insideNum`: Number of refined Ritz values within the interval, before removing copies.
    - `sinSearch`: Sine of angles between the desired eigenspace and the search subspace spanned by `S` before truncation at each iteration.
    - `sinRetain`: Sine of angles between the desired eigenspace and the approximate eigenspace spanned by the retained refined Ritz vectors.


With `maxItr=0`, no iteration is performed and `info.itrInfo` is absent.


## Related Publication

This implementation is based on the research described in: "Stable and effective restarting of the SS--RR and SS--RRR methods for computing partial eigenpairs of Hermitian matrices" by Zhongxiao Jia and Tianhang Liu.
