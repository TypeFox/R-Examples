library(easyVerification)
context('Error scores and skill scores')

fcst <- array(1.5, c(3,10))
obs <- rep(1, 3)

test_that('Mean error (bias)', {
  expect_is(EnsMe(fcst, obs), 'numeric')  
  expect_equal(EnsMe(fcst, obs), 0.5)
  expect_equal(EnsMess(fcst, fcst, obs), 0)
  expect_equal(EnsMess(fcst, fcst, obs), EnsErrorss(fcst, fcst, obs, type='me'))
})

test_that('Mean absolute error', {
  expect_is(EnsMe(fcst, obs), 'numeric')  
  expect_equal(EnsMae(fcst, obs), 0.5)
  expect_equal(EnsMaess(fcst, fcst, obs), 0)
})

test_that('Mean squared error', {
  expect_is(EnsMse(fcst, obs), 'numeric')  
  expect_equal(EnsMse(fcst, obs), 0.5**2)
  expect_equal(EnsMsess(fcst, fcst, obs), 0)
  expect_true(EnsMsess(fcst, fcst, obs) < EnsMsess(fcst, array(0, dim(fcst)), obs))
})

test_that('Root mean squared error', {
  expect_is(EnsRmse(fcst, obs), 'numeric')  
  expect_equal(EnsRmse(fcst, obs), 0.5)
  expect_equal(EnsRmse(fcst, obs), sqrt(EnsMse(fcst, obs)))
  expect_equal(EnsRmsess(fcst, fcst, obs), 0)
  expect_true(EnsRmsess(fcst, array(0, dim(fcst)), obs) < EnsMsess(fcst, array(0, dim(fcst)), obs))
})

