# TSCure

**Two-Step Estimation for Case-Cohort Studies with a Cure Fraction**

`TSCure` is an R package for semiparametric transformation mixture cure models under generalized case-cohort sampling. It implements inverse-probability-weighted (IPW) sieve estimation using an EM algorithm and an auxiliary-variable working-model update intended to improve estimation efficiency. It also provides multiplier-bootstrap inference, AIC-based model selection, and simulation tools.

The incidence component models the probability of being **uncured**. The transformation family includes proportional hazards (`r = 0`) and proportional odds (`r = 1`) models.

## Installation

`TSCure` is currently distributed as a source package. From a local package directory:

```r
install.packages("remotes")
remotes::install_local("/path/to/TSCure", dependencies = TRUE)
library(TSCure)
```

After the GitHub repository is published, a GitHub installation command can be added here.

## Quick start: simulate data and fit a model

```r
library(TSCure)

sim <- tscure_simulate(
  n = 500,
  r = 0,
  tau = 4.5,
  q1 = 0.2,
  q2 = 0.5,
  error_sd = 0.8,
  seed = 527
)

fit <- tscure_fit(
  data = sim,
  expensive = "U",
  auxiliary = "U_a",
  covariates = "Z",
  ind1 = "Ind1",
  ind2 = "Ind2",
  q1 = 0.2,
  q2 = 0.5,
  r = 0,
  n_knots = 1,
  degree = 1,
  nboot = 100,
  update = TRUE,
  seed = 527
)

coef(fit)
coef(fit, updated = TRUE)
summary(fit)
summary(fit, updated = TRUE)
vcov(fit, updated = TRUE)
fit$convergence
```

The bootstrap can be computationally intensive. For a quick smoke test, use `nboot = 0, update = FALSE`; this skips inference and the updated estimator. With `update = TRUE`, the covariance adjustment requires enough successful bootstrap fits to estimate an invertible working-model difference covariance matrix.

## Main functions

| Function | Purpose |
| --- | --- |
| `tscure_simulate()` | Generate a generalized case-cohort mixture cure dataset |
| `tscure_weights()` | Calculate inverse-probability sampling weights |
| `tscure_fit()` | Fit original IPW and optionally updated estimators |
| `tscure_select()` | Compare transformation parameters and spline knot counts by AIC |
| `tscure_study()` | Run repeated simulations; the current development version supports progress reporting and checkpoints |

## Proportional odds model

Set `r = 1` in both data generation and model fitting. For example:

```r
po_data <- tscure_simulate(n = 500, r = 1, tau = 11, seed = 527)
po_fit <- tscure_fit(
  po_data,
  expensive = "U", auxiliary = "U_a", covariates = "Z",
  ind1 = "Ind1", ind2 = "Ind2", q1 = 0.2, q2 = 0.5,
  r = 1, nboot = 0, update = FALSE
)
coef(po_fit)
```

## AIC model selection

```r
selection <- tscure_select(
  data = sim,
  r_grid = c(0, 0.5, 1),
  knot_grid = 1:3,
  expensive = "U", auxiliary = "U_a", covariates = "Z",
  ind1 = "Ind1", ind2 = "Ind2", q1 = 0.2, q2 = 0.5
)
selection$table
selection$selected
```

## Working with your own data

`tscure_fit()` accepts user-supplied column names for observed follow-up time (`time`), event indicator (`status`), expensive covariate (`expensive`), optional full-cohort auxiliary variable (`auxiliary`), and inexpensive covariates (`covariates`). Sampling indicators (`ind1`, `ind2`) and probabilities (`q1`, `q2`) can be supplied for generalized case-cohort designs. Alternatively, provide an observation-weight vector using `weights`. If neither sampling indicators nor weights are provided, unit weights are used.

## Simulation validation

Preliminary checks of the development implementation used representative manuscript settings with 100 simulated cohorts and 100 bootstrap replications per cohort, at `n = 2000`. The PH scenario used `r = 0, tau = 4.5`, and the PO scenario used `r = 1, tau = 11`, with `q1 = 0.2`, `q2 = 0.5`, and auxiliary-variable noise SD `0.8`. The updated estimator showed lower empirical variance than the original IPW estimator for all five reported regression parameters in both scenarios. These checks are **not** exact reproductions of the publication's tables, which used more simulation and bootstrap replications. The PH experiment also showed below-nominal coverage for some updated latency coefficients, warranting caution when interpreting finite-sample inference.

The development package has passed `R CMD check` with zero errors, warnings, and notes in the maintainer's local R environment. Recheck after any further source or documentation changes.

## Reference

Qingning Zhou and Xu Cao. *Improving Estimation Efficiency for Case-cohort Studies with a Cure Fraction*. Biometrics. Please consult the published article for the definitive bibliographic details and theoretical results.

## Development status

Version `0.1.0` is an initial research-software release candidate. Users should inspect convergence diagnostics and assess performance under their own sampling and censoring settings. Exact numerical agreement with the published simulation tables is not claimed.

## License

MIT (see `LICENSE` in the source package).
