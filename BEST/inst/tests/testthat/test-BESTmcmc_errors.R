
# Tests for BESTmcmc

context("BESTmcmc_errors")

test_that("BESTmcmc gives sensible errors",  {
  expect_error(BESTmcmc(c(1:4, NA), 1:4, numSavedSteps = 9, burnInSteps = 1),
    "The input data include NA or Inf")
  expect_error(BESTmcmc(1:4, 1:4, 9),
    "'priors' is now the 3rd argument")
  expect_error(BESTmcmc(1:4, 1:4, 9),
    "it must be a list")
  expect_error(BESTmcmc(1:4, 1:4, priors="A"),
    "'priors' must be a list")
  expect_error(BESTmcmc(1:4, 1:4, list(nonsense=8)),
    "Invalid items in prior specification")
  expect_error(BESTmcmc(1:4, 1:4, list(muSD=-1)),
    "muSD must be > 0")
  expect_error(BESTmcmc(1:4, 1:4, list(sigmaMode=-1)),
    "gamma prior must be > 0")
  expect_error(BESTmcmc(1:4, 1:4, list(nuSD=-1)),
    "gamma prior must be > 0")
  expect_error(BESTmcmc(4, 4, list()),
    "data must include at least 2")
  expect_error(BESTmcmc(4, 4),
    "data must include at least 2")
})
