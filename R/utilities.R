
# ============================================================
# TSCure: Simulation and Utility Functions
# Version 0.1.0
# ============================================================


#' Simulate a Generalized Case-Cohort Mixture Cure Dataset
#'
#' Generate survival data under a semiparametric transformation
#' mixture cure model with generalized case-cohort sampling.
#'
#' @param n Cohort size.
#' @param alpha_A,gamma_A Latency regression coefficients.
#' @param alpha_I,gamma_I,lambda_I Incidence regression
#'   coefficients and intercept.
#' @param tau Administrative study duration.
#' @param r Transformation parameter (0 = PH, 1 = PO).
#' @param q1,q2 Two-phase sampling probabilities.
#' @param error_sd Standard deviation of auxiliary-variable noise.
#' @param seed Optional random seed.
#'
#' @return A data frame containing A, T, C, Y, D, U, Z,
#'   Ind1, Ind2, and U_a.
#' @export
tscure_simulate <- function(
    n = 2000,
    alpha_A = log(2),
    gamma_A = -0.5,
    alpha_I = log(2),
    gamma_I = -0.5,
    lambda_I = -0.5,
    tau = 4.5,
    r = 0,
    q1 = 0.2,
    q2 = 0.5,
    error_sd = 0.8,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  stopifnot(
    n > 1,
    tau > 0,
    r >= 0,
    error_sd >= 0
  )
  
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  
  U <- z1
  Z <- 0.5 * z1 + sqrt(0.75) * z2
  
  U_a <- U + stats::rnorm(n, sd = error_sd)
  
  pi <- stats::plogis(
    lambda_I + alpha_I * U + gamma_I * Z
  )
  
  A <- stats::rbinom(n, 1, pi)
  
  ru <- stats::runif(n)
  h <- -log1p(-ru)
  
  target <- if (r == 0) {
    h
  } else {
    expm1(r * h) / r
  }
  
  # Baseline cumulative hazard: (t^2 + t) / 2
  v <- target / exp(alpha_A * U + gamma_A * Z)
  
  T <- ifelse(
    A == 1,
    (-1 + sqrt(1 + 8 * v)) / 2,
    Inf
  )
  
  C <- pmin(
    stats::runif(n, 0, 5 * tau / 4),
    tau
  )
  
  Y <- pmin(T, C)
  D <- as.integer(T <= C)
  
  Ind1 <- stats::rbinom(n, 1, q1)
  
  Ind2 <- integer(n)
  
  eligible <- which(D == 1 & Ind1 == 0)
  
  Ind2[eligible] <- stats::rbinom(
    length(eligible),
    1,
    q2
  )
  
  data.frame(
    A, T, C, Y, D, U, Z,
    Ind1, Ind2, U_a
  )
}


#' Select Transformation Parameter and Spline Knots by AIC
#'
#' Fit candidate mixture cure models using different
#' transformation parameters and numbers of interior knots.
#'
#' @param data Dataset.
#' @param r_grid Candidate transformation parameters.
#' @param knot_grid Candidate interior-knot counts.
#' @param ... Additional arguments passed to tscure_fit,
#'   excluding r, n_knots, nboot, and update.
#'
#' @return A list containing the AIC table, selected model,
#'   and best fitted model.
#' @export
tscure_select <- function(
    data,
    r_grid = seq(0, 3, 0.5),
    knot_grid = 1:4,
    ...
) {
  
  grid <- expand.grid(
    r = r_grid,
    n_knots = knot_grid
  )
  
  fits <- vector("list", nrow(grid))
  
  grid$AIC <- NA_real_
  grid$converged <- FALSE
  
  for (i in seq_len(nrow(grid))) {
    
    fit <- tryCatch(
      tscure_fit(
        data,
        ...,
        r = grid$r[i],
        n_knots = grid$n_knots[i],
        nboot = 0,
        update = FALSE
      ),
      error = function(e) NULL
    )
    
    fits[[i]] <- fit
    
    if (!is.null(fit)) {
      
      grid$AIC[i] <-
        -2 * fit$logLik +
        2 * length(fit$original_full)
      
      grid$converged[i] <-
        isTRUE(fit$convergence$original)
    }
  }
  
  ok <- which(
    grid$converged &
      is.finite(grid$AIC)
  )
  
  if (!length(ok)) {
    stop("No candidate models converged.")
  }
  
  best <- ok[which.min(grid$AIC[ok])]
  
  list(
    table = grid[order(grid$AIC), , drop = FALSE],
    best = fits[[best]],
    selected = grid[best, , drop = FALSE]
  )
}


#' Run a Simulation Study Comparing IPW and Updated Estimators
#'
#' Repeatedly generate generalized case-cohort datasets and
#' estimate mixture cure model parameters using the original
#' IPW estimator and the updated estimator.
#'
#' Supports progress reporting and automatic checkpoint saving.
#'
#' @param nsim Number of simulated cohorts.
#' @param n Cohort size.
#' @param nboot Number of bootstrap replications per cohort.
#' @param seed Random seed.
#' @param r Transformation parameter (0 = PH, 1 = PO).
#' @param tau Administrative study duration.
#' @param q1 First-phase sampling probability.
#' @param q2 Second-phase case sampling probability.
#' @param error_sd Auxiliary-variable error standard deviation.
#' @param n_knots Number of interior spline knots.
#' @param degree Spline degree.
#' @param max_iter Maximum number of EM iterations.
#' @param optim_maxit Maximum optimizer iterations.
#' @param checkpoint Optional path for saving intermediate results.
#' @param progress Whether to display progress messages.
#'
#' @return A list containing fitted models, error messages,
#'   original estimates, updated estimates, number of completed
#'   simulations, and elapsed computation time.
#' @export
tscure_study <- function(
    nsim = 10,
    n = 2000,
    nboot = 5,
    seed = 5,
    r = 0,
    tau = 4.5,
    q1 = 0.2,
    q2 = 0.5,
    error_sd = 0.8,
    n_knots = 1,
    degree = 1,
    max_iter = 5000,
    optim_maxit = 5000,
    checkpoint = NULL,
    progress = TRUE
) {
  
  stopifnot(
    length(nsim) == 1L,
    is.finite(nsim),
    nsim >= 1,
    nsim == as.integer(nsim)
  )
  
  set.seed(seed)
  
  fits <- vector("list", nsim)
  errors <- rep(NA_character_, nsim)
  
  start_time <- Sys.time()
  
  # Construct results from completed simulations
  make_results <- function(k) {
    
    current <- fits[seq_len(k)]
    
    good <- Filter(
      Negate(is.null),
      current
    )
    
    original <- if (length(good)) {
      do.call(
        rbind,
        lapply(good, `[[`, "original")
      )
    } else {
      NULL
    }
    
    updated <- if (
      length(good) &&
      all(vapply(
        good,
        function(x) !is.null(x$updated),
        logical(1)
      ))
    ) {
      do.call(
        rbind,
        lapply(good, `[[`, "updated")
      )
    } else {
      NULL
    }
    
    list(
      fits = current,
      errors = errors[seq_len(k)],
      original = original,
      updated = updated,
      completed = k,
      nsim = nsim,
      elapsed_seconds = as.numeric(
        difftime(
          Sys.time(),
          start_time,
          units = "secs"
        )
      )
    )
  }
  
  # Main simulation loop
  for (k in seq_len(nsim)) {
    
    d <- tscure_simulate(
      n = n,
      r = r,
      tau = tau,
      q1 = q1,
      q2 = q2,
      error_sd = error_sd
    )
    
    fits[[k]] <- tryCatch(
      tscure_fit(
        d,
        ind1 = "Ind1",
        ind2 = "Ind2",
        q1 = q1,
        q2 = q2,
        r = r,
        n_knots = n_knots,
        degree = degree,
        nboot = nboot,
        max_iter = max_iter,
        optim_maxit = optim_maxit
      ),
      error = function(e) {
        errors[k] <<- conditionMessage(e)
        NULL
      }
    )
    
    result <- make_results(k)
    
    # Save checkpoint after every completed simulation
    if (!is.null(checkpoint)) {
      
      tmp <- tempfile(
        pattern = ".tscure_",
        tmpdir = dirname(checkpoint),
        fileext = ".rds"
      )
      
      saveRDS(result, tmp)
      
      if (!file.rename(tmp, checkpoint)) {
        if (file.exists(checkpoint)) {
          unlink(checkpoint)
        }
        
        if (!file.rename(tmp, checkpoint)) {
          stop(
            "Could not save checkpoint: ",
            checkpoint
          )
        }
      }
    }
    
    # Display progress
    if (isTRUE(progress)) {
      
      elapsed <- result$elapsed_seconds
      
      cat(
        sprintf(
          "Simulation %d/%d completed | %.1f minutes elapsed\n",
          k,
          nsim,
          elapsed / 60
        )
      )
      
      utils::flush.console()
    }
  }
  
  make_results(nsim)
}
