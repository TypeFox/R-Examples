context("Testing moment estimators \n")
set.seed(40)
nobs <- 1e3

yy <- rcauchy(n = 1000)

test_that("skewness is approximately zero for a Normal distribution", {
  skew.sim <- replicate(1000, skewness(rnorm(nobs)))
  expect_equal(mean(skew.sim), 0, tol = 0.01)
})

test_that("kurtosis is approximately 3 for a Normal distribution", {
  kurt.sim <- replicate(1000, kurtosis(rnorm(nobs)))
  expect_equal(mean(kurt.sim), 3, tol = 0.01)
})


test_that("kurtosis is computing the fourth central moment", {
  expect_equal(kurtosis(yy), mean((yy - mean(yy))^4) / var(yy)^2, tol = 0.01)
})

test_that("skewness is computing the fourth central moment", {
  expect_equal(skewness(yy), mean((yy - mean(yy))^3) / var(yy)^1.5, tol = 0.01)
})
