context("Testing arima_pi")

set.seed(123)
x <- rnorm(10)

test_that("bogus arguments throw error",{
  expect_error(arima_pi(x, phi = -2))
  expect_error(arima_pi(x, order = "ar"))
  expect_error(arima_pi(x, xreg = 1))
  expect_error(arima_pi(x, n_ahead = -1))
  expect_error(arima_pi(x, level = 2))
  expect_error(arima_pi(x, median = "true"))
  expect_error(arima_pi(x, prior = "custom"))
  expect_error(arima_pi(x, prior = "custom", custom_prior = "f"))
  expect_error(arima_pi(x, nsim = 0))
})

test_that("output of arima_pi is of correct size and form",{
  y <- ts(rnorm(10), start = 2000, frequency = 12)
  pred <- arima_pi(y, order = c(1, 0, 0), n_ahead = 5, nsim = 50)
  expect_identical(dim(pred), c(5L, 5L))
  expect_identical(class(pred), c("mts", "ts", "matrix"))
  expect_identical(frequency(pred), frequency(y))
  expect_identical(start(pred), end(y)+c(0, 1))
})

test_that("same seeds give same results",{
    set.seed(1)
    pred1 <- arima_pi(x, c(1, 0, 0), nsim = 50)
    set.seed(1)
    pred2 <- arima_pi(x, c(1, 0, 0), nsim = 50)
    expect_identical(pred1, pred2)
})

test_that("larger nsim gives smaller se",{
  set.seed(1)
  pred1 <- arima_pi(x, c(1, 0, 0), nsim = 50)
  set.seed(1)
  pred2 <- arima_pi(x, c(1, 0, 0), nsim = 100)
  expect_gt(as.numeric(pred1[, "se_upr"]), as.numeric(pred2[, "se_upr"]))
})

test_that("arima_pi with uniform prior gives same results each time",{
  set.seed(1)
  pred <- arima_pi(lh, c(1, 0, 0), nsim = 50)
  expect_equal(pred[1,"median"], 2.707872101, tol = 1e-5, check.attributes = FALSE)
  expect_equal(pred[1,"lwr"], 1.809644512, tol = 1e-5, check.attributes = FALSE)
  expect_equal(pred[1,"upr"], 3.606841626, tol = 1e-5, check.attributes = FALSE)
  expect_equal(pred[1,"se_lwr"], 0.01355273, tol = 1e-5, check.attributes = FALSE)
  expect_equal(pred[1,"se_upr"], 0.0157178, tol = 1e-5, check.attributes = FALSE)
})
test_that("arima_pi with jeffreys gives same results each time",{
  set.seed(1)
  pred <- arima_pi(lh, c(1, 0, 0), nsim = 50, prior = "approx_marginal")
  expect_equal(pred[1,"lwr"], 1.813243, tol = 1e-5, check.attributes = FALSE)
  pred <- arima_pi(lh, c(1, 0, 0), nsim = 50, prior = "approx_joint")
  expect_equal(pred[1,"lwr"], 1.759795, tol = 1e-5, check.attributes = FALSE)
  pred <- arima_pi(lh, c(1, 0, 0), nsim = 50, prior = "exact_marginal")
  expect_equal(pred[1,"lwr"], 1.795299, tol = 1e-5, check.attributes = FALSE)
  pred <- arima_pi(lh, c(1, 0, 0), nsim = 50, prior = "exact_joint")
  expect_equal(pred[1,"lwr"], 1.7381, tol = 1e-5, check.attributes = FALSE)
})
