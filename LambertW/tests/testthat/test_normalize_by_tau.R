context("Testing normalize by tau function\n")

random.data <- rnorm(100)

wrong.taus <- list(c("mu" = 0, "sigma" = 1),
                   c("mean" = 0, "sd" = 1),
                   c(0, 1),
                   c("mu_x" = 0))

correct.taus <- list("identity" = c("mu_x" = 0, "sigma_x" = 1),
                     "identity2" = c("mu_x" = 0, "sigma_x" = 1, "gamma" = 1),
                     "times2" = c("mu_x" = 0, "sigma_x" = 2),
                     "plus1times2" = c("mu_x" = -1, "sigma_x" = 2))

test_that("throw errors for incorrect input arguments", {
  
  for (tt in correct.taus) {
    expect_true(is.numeric(normalize_by_tau(random.data, tt)))
    # default is inverse = FALSE
    expect_equal(normalize_by_tau(random.data, tt),
                 normalize_by_tau(random.data, tt, inverse = FALSE))
  }
  for (tt in wrong.taus) {
    expect_error(normalize_by_tau(random.data, tt))
  }
  expect_error(normalize_by_tau())
  expect_error(normalize_by_tau(random.data, correct.tau, 0))
})

test_that("it computes the correct transformation", {
  
  for (nn in names(correct.taus)) {
    tau.tmp <- correct.taus[[nn]]
    for (ll in c(FALSE, TRUE)) {
      if (grepl("identity", nn)) {
        expect_equal(normalize_by_tau(random.data, tau.tmp, inverse = ll),
                     random.data)
      }
    }
    
    expect_equal(normalize_by_tau(random.data, tau.tmp, inverse = FALSE),
                 (random.data - tau.tmp["mu_x"]) / tau.tmp["sigma_x"])
    # inverse is actually the inverse
    norm.data <- normalize_by_tau(random.data, tau.tmp, inverse = FALSE)
    expect_equal(normalize_by_tau(norm.data, tau.tmp, inverse = TRUE),
                 random.data)
  }
})
