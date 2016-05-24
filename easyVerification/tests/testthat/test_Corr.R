library(easyVerification)
context('Test correlation')

obs <- rnorm(10)
fcst <- array(2, c(10,3)) + obs
fcst2 <- array(rnorm(30), c(10,3)) + obs

test_that('output type and results', {
  expect_is(EnsCorr(fcst, obs), 'numeric')
  expect_equal(EnsCorr(fcst, obs), 1)
  expect_equal(EnsCorr(fcst, obs), cor(rowMeans(fcst), obs))
  expect_true(EnsCorr(fcst2 + obs, obs) >= EnsCorr(fcst2, obs))
})