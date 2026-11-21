#' Generalized case-cohort inverse probability weights
#' @param status Binary event indicator.
#' @param ind1 Subcohort inclusion indicator.
#' @param ind2 Additional case inclusion indicator.
#' @param q1 Subcohort sampling probability.
#' @param q2 Additional case sampling probability outside subcohort.
#' @return Numeric weights (zero for unselected subjects).
#' @export
#' @examples
#' tscure_weights(c(0, 1, 1), c(1, 0, 0), c(0, 1, 0), .2, .5)
tscure_weights <- function(status, ind1, ind2, q1, q2) {
  if (length(status) != length(ind1) || length(status) != length(ind2)) stop('Length mismatch.')
  if (anyNA(c(status, ind1, ind2)) ||
      !all(status %in% 0:1) || !all(ind1 %in% 0:1) || !all(ind2 %in% 0:1))
    stop('Indicators must be binary without missing values.')
  if (!is.finite(q1) || !is.finite(q2) || q1 <= 0 || q1 > 1 || q2 < 0 || q2 > 1)
    stop('Invalid sampling probabilities.')
  if (any(ind2 == 1 & (status != 1 | ind1 == 1))) stop('ind2 must select only cases outside the subcohort.')
  selected <- (ind1 == 1 | ind2 == 1)
  p <- ifelse(status == 1, q1 + (1 - q1) * q2, q1)
  as.numeric(selected) / p
}
#' Fit a semiparametric transformation mixture cure model
#'
#' Uses an inverse-probability weighted EM algorithm, and optionally a
#' working model on the full cohort to update the regression coefficients.
#' The incidence component models the probability of being uncured.
#'
#' @param data Data frame containing time, event and covariates.
#' @param time,status Column names for observed time and event (0/1).
#' @param expensive Column name of expensive covariate.
#' @param auxiliary Optional column name of full-cohort proxy. If NULL, the
#'   working model uses only `covariates`.
#' @param covariates Character vector of inexpensive adjustment covariates.
#' @param ind1,ind2 Optional sampling-indicator column names.
#' @param q1,q2 Sampling probabilities; required when indicators are supplied.
#' @param weights Optional numeric observation weights; overrides sampling indicators.
#' @param r Nonnegative logarithmic transformation parameter (0 = PH, 1 = PO).
#' @param n_knots Number of interior I-spline knots.
#' @param degree Spline degree (1 reproduces manuscript settings).
#' @param nboot Number of multiplier bootstrap replicates; 0 skips inference.
#' @param update Whether to fit the working models and update the estimator.
#' @param seed Optional random seed for multiplier weights.
#' @param max_iter,tol,optim_maxit EM iterations, convergence tolerance and
#'   maximum Nelder-Mead iterations per M-step.
#' @return Object of class `tscure_fit` with original and updated estimates,
#'   variance matrices and convergence diagnostics.
#' @export
#' @examples
#' d <- tscure_simulate(n = 120, seed = 10)
#' fit <- tscure_fit(d, expensive = 'U', auxiliary = 'U_a',
#'                   covariates = 'Z', ind1 = 'Ind1', ind2 = 'Ind2',
#'                   q1 = .2, q2 = .5, nboot = 0, update = FALSE,
#'                   max_iter = 10, optim_maxit = 50)
tscure_fit <- function(data, time = 'Y', status = 'D', expensive = 'U',
                       auxiliary = 'U_a', covariates = 'Z', ind1 = NULL,
                       ind2 = NULL, q1 = NULL, q2 = NULL, weights = NULL,
                       r = 0, n_knots = 1, degree = 1, nboot = 0,
                       update = TRUE, seed = NULL, max_iter = 5000,
                       tol = 1e-4, optim_maxit = 5000) {
  stopifnot(is.data.frame(data), is.character(covariates), length(expensive) == 1,
            length(r) == 1, is.finite(r), r >= 0, n_knots >= 0,
            n_knots == as.integer(n_knots), degree >= 1, degree == as.integer(degree),
            nboot >= 0, nboot == as.integer(nboot), max_iter >= 1, tol > 0)
  needed <- unique(c(time, status, expensive, covariates, if (update) auxiliary,
                     ind1, ind2))
  if (!all(needed %in% names(data))) stop('Missing data columns: ', paste(setdiff(needed, names(data)), collapse = ', '))
  n <- nrow(data)
  if (is.null(weights)) {
    if (is.null(ind1) != is.null(ind2)) stop('Provide both ind1 and ind2.')
    if (is.null(ind1)) weights <- rep(1, n) else {
      if (is.null(q1) || is.null(q2)) stop('Provide q1 and q2.')
      weights <- tscure_weights(data[[status]], data[[ind1]], data[[ind2]], q1, q2)
    }
  }
  if (length(weights) != n || anyNA(weights) || any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0))
    stop('Invalid weights.')
  y <- data[[time]]
  if (anyNA(y) || any(!is.finite(y)) || any(y < 0)) stop('Invalid times.')
  if (anyNA(data[[status]]) || !all(data[[status]] %in% 0:1)) stop('Invalid status.')
  if (!any(data[[status]] == 1 & weights > 0)) stop('No weighted observed events.')
  knots <- if (n_knots > 0) as.numeric(stats::quantile(y, probs = seq_len(n_knots) / (n_knots + 1))) else numeric()
  boundary <- c(0, max(y))
  active <- which(weights > 0)
  d_orig <- .tscure_design(data[active, , drop = FALSE], time, status,
                           c(expensive, covariates), knots, degree, boundary)
  original <- .tscure_em(d_orig, weights[active], r, max_iter, tol, optim_maxit)
  work_cc <- work_full <- NULL
  if (update) {
    work_names <- c(if (!is.null(auxiliary)) auxiliary, covariates)
    if (!length(work_names)) stop('Working model requires auxiliary or covariates.')
    d_cc <- .tscure_design(data[active, , drop = FALSE], time, status,
                            work_names, knots, degree, boundary)
    d_full <- .tscure_design(data, time, status, work_names, knots, degree, boundary)
    work_cc <- .tscure_em(d_cc, weights[active], r, max_iter, tol, optim_maxit)
    work_full <- .tscure_em(d_full, rep(1, n), r, max_iter, tol, optim_maxit)
  }
  p <- ncol(d_orig$x); m <- ncol(d_orig$basis)
  idx <- c(seq_len(p), p + m + seq_len(p + 1))
  param_names <- c(paste0('latency:', colnames(d_orig$x)), 'incidence:(Intercept)',
                   paste0('incidence:', colnames(d_orig$x)))
  # Same order as internal parameter vector after dropping spline coefficients.
  point <- original$par[idx]; names(point) <- param_names
  out <- list(original = point, updated = NULL, original_full = original$par,
              working_casecohort = work_cc, working_full = work_full,
              variance_original = NULL, variance_updated = NULL,
              bootstrap = NULL, knots = knots, degree = degree, r = r,
              convergence = list(original = original$converged,
                working_casecohort = if (update) work_cc$converged else NULL,
                working_full = if (update) work_full$converged else NULL),
              logLik = -original$negloglik, n = n, n_selected = length(active),
              nboot = nboot, call = match.call())
  if (update) {
    # Working-model regression coefficients can have a different dimension.
    pw <- ncol(d_cc$x); idxw <- c(seq_len(pw), pw + m + seq_len(pw + 1))
    diff <- work_cc$par[idxw] - work_full$par[idxw]
  }
  if (nboot > 0) {
    if (nboot < 2) stop('At least 2 bootstrap replicates required.')
    if (!is.null(seed)) {if (exists('.Random.seed', .GlobalEnv)) old_seed <- get('.Random.seed', .GlobalEnv) else old_seed <- NULL
      on.exit({if (is.null(old_seed)) {if (exists('.Random.seed', .GlobalEnv)) rm('.Random.seed', envir = .GlobalEnv)} else assign('.Random.seed', old_seed, envir = .GlobalEnv)}, add = TRUE)
      set.seed(seed)}
    bp <- matrix(NA_real_, nboot, length(idx))
    if (update) bd <- matrix(NA_real_, nboot, length(diff))
    boot_ok <- rep(FALSE, nboot)
    for (b in seq_len(nboot)) {
      u <- stats::rexp(n)
      a <- .tscure_em(d_orig, weights[active] * u[active], r, max_iter, tol, optim_maxit)
      if (update) {
        wc <- .tscure_em(d_cc, weights[active] * u[active], r, max_iter, tol, optim_maxit)
        wf <- .tscure_em(d_full, u, r, max_iter, tol, optim_maxit)
        boot_ok[b] <- a$converged && wc$converged && wf$converged
        bd[b, ] <- wc$par[idxw] - wf$par[idxw]
      } else boot_ok[b] <- a$converged
      bp[b, ] <- a$par[idx]
    }
    if (!all(boot_ok)) warning(sum(!boot_ok), ' bootstrap fits did not meet convergence tolerance; inspect diagnostics.')
    good <- boot_ok & stats::complete.cases(bp) & if (update) stats::complete.cases(bd) else TRUE
    if (sum(good) < 2) stop('Too few converged bootstrap fits to estimate variance.')
    bp <- bp[good, , drop = FALSE]
    out$variance_original <- stats::cov(bp)
    dimnames(out$variance_original) <- list(param_names, param_names)
    if (update) {
      bd <- bd[good, , drop = FALSE]
      s12 <- stats::cov(bp, bd); s22 <- stats::cov(bd)
      if (rcond(s22) < 1e-10) stop('Working-model difference covariance is singular/ill-conditioned; increase nboot or revise working model.')
      adj <- s12 %*% solve(s22)
      out$updated <- drop(point - adj %*% diff)
      names(out$updated) <- param_names
      out$variance_updated <- out$variance_original - adj %*% t(s12)
      dimnames(out$variance_updated) <- list(param_names, param_names)
    }
    out$bootstrap <- list(original = bp, difference = if (update) bd else NULL,
                          converged = boot_ok, n_success = sum(good))
  }
  class(out) <- 'tscure_fit'
  out
}
#' @export
coef.tscure_fit <- function(object, updated = FALSE, ...) {
  if (updated) {if (is.null(object$updated)) stop('Updated coefficients unavailable.'); return(object$updated)}
  object$original
}
#' @export
vcov.tscure_fit <- function(object, updated = FALSE, ...) {
  v <- if (updated) object$variance_updated else object$variance_original
  if (is.null(v)) stop('Bootstrap covariance unavailable; fit with nboot >= 2.')
  v
}
#' @export
summary.tscure_fit <- function(object, updated = FALSE, ...) {
  b <- coef(object, updated = updated)
  v <- if (updated) object$variance_updated else object$variance_original
  se <- if (is.null(v)) rep(NA_real_, length(b)) else sqrt(pmax(diag(v), 0))
  tab <- cbind(Estimate = b, `Std. Error` = se,
               `z value` = b / se, `Pr(>|z|)` = 2 * stats::pnorm(-abs(b / se)))
  out <- list(coefficients = tab, updated = updated, convergence = object$convergence,
              n = object$n, n_selected = object$n_selected)
  class(out) <- 'summary.tscure_fit'; out
}
#' @export
print.summary.tscure_fit <- function(x, ...) {print(x$coefficients); invisible(x)}
#' @export
print.tscure_fit <- function(x, ...) {
  cat('TSCure fit: n =', x$n, ', selected =', x$n_selected,
      ', r =', x$r, ', knots =', length(x$knots), '\n')
  print(x$original)
  if (!is.null(x$updated)) {cat('Updated estimates:\n'); print(x$updated)}
  invisible(x)
}
