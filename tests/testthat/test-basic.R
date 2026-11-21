
# Tests for TSCure

test_that("Simulation generates valid data", {
  
  d <- tscure_simulate(
    n = 100,
    seed = 527
  )
  
  expect_equal(nrow(d), 100)
  
  expect_true(all(c(
    "A", "T", "C", "Y", "D",
    "U", "Z", "Ind1", "Ind2", "U_a"
  ) %in% names(d)))
  
  expect_true(all(d$D %in% c(0, 1)))
  expect_true(all(d$Ind1 %in% c(0, 1)))
  expect_true(all(d$Ind2 %in% c(0, 1)))
  
  expect_equal(d$Y, pmin(d$T, d$C))
  
  expect_true(all(
    d$Ind2 == 0 | (d$D == 1 & d$Ind1 == 0)
  ))
})


test_that("Simulation is reproducible", {
  
  d1 <- tscure_simulate(
    n = 100,
    seed = 527
  )
  
  d2 <- tscure_simulate(
    n = 100,
    seed = 527
  )
  
  expect_equal(d1, d2)
})


test_that("Sampling weights are valid", {
  
  d <- tscure_simulate(
    n = 200,
    seed = 527
  )
  
  w <- tscure_weights(
    d$D,
    d$Ind1,
    d$Ind2,
    q1 = 0.2,
    q2 = 0.5
  )
  
  expect_length(w, nrow(d))
  expect_true(all(is.finite(w)))
  expect_true(all(w >= 0))
  
  expect_equal(
    w[d$D == 0 & d$Ind1 == 1],
    rep(5, sum(d$D == 0 & d$Ind1 == 1))
  )
})


test_that("PH model fitting converges", {
  
  d <- tscure_simulate(
    n = 300,
    r = 0,
    seed = 527
  )
  
  fit <- tscure_fit(
    data = d,
    ind1 = "Ind1",
    ind2 = "Ind2",
    q1 = 0.2,
    q2 = 0.5,
    r = 0,
    nboot = 0,
    update = FALSE
  )
  
  expect_true(
    isTRUE(fit$convergence$original)
  )
  
  expect_true(
    all(is.finite(fit$original))
  )
  
  expect_length(fit$original, 5)
})


test_that("PO model fitting returns finite estimates", {
  
  d <- tscure_simulate(
    n = 300,
    r = 1,
    seed = 527
  )
  
  fit <- tscure_fit(
    data = d,
    ind1 = "Ind1",
    ind2 = "Ind2",
    q1 = 0.2,
    q2 = 0.5,
    r = 1,
    nboot = 0,
    update = FALSE
  )
  
  expect_true(
    all(is.finite(fit$original))
  )
  
  expect_length(fit$original, 5)
})
