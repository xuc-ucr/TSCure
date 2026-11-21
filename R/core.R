# Internal model implementation. Incidence is probability of being UNCURED.
#' @importFrom stats coef plogis qlogis
NULL

.tscure_transform <- function(x, r) {
  if (r == 0) return(x)
  log1p(r * x) / r
}
.tscure_log_deriv <- function(x, r) {
  if (r == 0) return(rep(0, length(x)))
  -log1p(r * x)
}
.tscure_split <- function(par, p, m) {
  list(latency = par[seq_len(p)],
       eta = par[p + seq_len(m)],
       incidence = par[p + m + seq_len(p + 1L)])
}
.tscure_design <- function(data, time, status, covariates, knots, degree, boundary) {
  y <- as.numeric(data[[time]])
  d <- as.numeric(data[[status]])
  if (anyNA(y) || anyNA(d) || any(y < 0) || !all(d %in% c(0, 1)))
    stop('Time must be nonnegative and status must be 0/1 with no missing values.')
  if (!any(d == 1)) stop('At least one observed event is required.')
  x <- as.matrix(data[, covariates, drop = FALSE])
  storage.mode(x) <- 'double'
  if (anyNA(x) || any(!is.finite(x))) stop('Covariates must be finite and observed for the fitted rows.')
  ib <- splines2::ibs(y, knots = knots, degree = degree,
                       intercept = TRUE, Boundary.knots = boundary)
  list(y = y, d = d, x = x, basis = unclass(ib),
       deriv = unclass(stats::deriv(ib)),
       max_event = max(y[d == 1]))
}
.tscure_components <- function(par, design, r) {
  p <- ncol(design$x); m <- ncol(design$basis)
  v <- .tscure_split(par, p, m)
  haz0 <- drop(design$deriv %*% exp(v$eta))
  cum0 <- drop(design$basis %*% exp(v$eta))
  lin_a <- drop(design$x %*% v$latency)
  lin_i <- v$incidence[1] + drop(design$x %*% v$incidence[-1])
  z <- cum0 * exp(lin_a)
  list(haz0 = haz0, cum0 = cum0, lin_a = lin_a,
       pi = plogis(lin_i), z = z, log_surv = -.tscure_transform(z, r),
       log_dg = .tscure_log_deriv(z, r))
}
.tscure_expect <- function(par, design, r) {
  a <- .tscure_components(par, design, r)
  # Taylor (1995) zero-tail convention used in the original manuscript code.
  ea <- design$d + (1 - design$d) * plogis(qlogis(a$pi) + a$log_surv)
  ea[design$d == 0 & design$y > design$max_event] <- 0
  ea
}
.tscure_q <- function(par, design, r, ea, weights) {
  a <- .tscure_components(par, design, r)
  if (any(!is.finite(a$haz0)) || any(a$haz0 <= 0) || any(!is.finite(a$z))) return(1e100)
  lp <- qlogis(a$pi)
  # Stable logistic log probabilities.
  log_inc <- ea * (-pmax(lp, 0) - log1p(exp(-abs(lp))) + lp) +
    (1 - ea) * (-pmax(lp, 0) - log1p(exp(-abs(lp))))
  log_lat <- design$d * (log(a$haz0) + a$lin_a + a$log_dg) + ea * a$log_surv
  out <- -sum(weights * (log_inc + log_lat))
  if (is.finite(out)) out else 1e100
}
.tscure_observed <- function(par, design, r, weights) {
  a <- .tscure_components(par, design, r)
  if (any(a$haz0 <= 0)) return(Inf)
  log_pi <- -pmax(-qlogis(a$pi), 0) - log1p(exp(-abs(qlogis(a$pi))))
  log_event <- log_pi + log(a$haz0) + a$lin_a + a$log_dg + a$log_surv
  log_cens <- log1p(-a$pi * (-expm1(a$log_surv)))
  out <- -sum(weights * ifelse(design$d == 1, log_event, log_cens))
  if (is.finite(out)) out else Inf
}
.tscure_em <- function(design, weights, r, max_iter, tol, optim_maxit, initial = NULL) {
  p <- ncol(design$x); m <- ncol(design$basis)
  if (is.null(initial)) initial <- rep(0, 2 * p + m + 1L)
  par <- initial; conv <- FALSE; opt_codes <- integer()
  for (it in seq_len(max_iter)) {
    ea <- .tscure_expect(par, design, r)
    opt <- stats::optim(par, .tscure_q, design = design, r = r,
                        ea = ea, weights = weights, method = 'Nelder-Mead',
                        control = list(maxit = optim_maxit))
    opt_codes <- c(opt_codes, opt$convergence)
    if (any(!is.finite(opt$par))) break
    delta <- sum(abs(opt$par - par))
    par <- opt$par
    if (delta <= tol && opt$convergence == 0) { conv <- TRUE; break }
  }
  list(par = par, converged = conv, iterations = length(opt_codes),
       optimizer_codes = opt_codes, negloglik = .tscure_observed(par, design, r, weights))
}
