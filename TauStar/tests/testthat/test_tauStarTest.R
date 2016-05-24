library(TauStar)
context("Testing the tauStarTest functionality.")

test_that("tauStarTest with continuous data.", {
  set.seed(238)
  sims = 5
  n = 100
  for (i in 1:sims) {
    x = rnorm(n)
    y = rnorm(n)
    tval = tStar(x, y)
    a = tauStarTest(x, y, mode = "continuous")
    expect_equal(a$x, x)
    expect_equal(a$y, y)
    expect_equal(a$tStar, tval)
    expect_equal(a$pVal, 1 - pHoeffInd(n * tval))
    expect_equal(a$mode, "continuous")
  }

  x = rnorm(n)
  y = rnorm(n)
  tval = tStar(x, y)
  a = tauStarTest(x, y)
  expect_equal(a$x, x)
  expect_equal(a$y, y)
  expect_equal(a$tStar, tval)
  expect_equal(a$pVal, 1 - pHoeffInd(n * tval))
  expect_equal(a$mode, "auto-continuous")

  expect_error(tauStarTest(rnorm(10), rnorm(9)))
})

test_that("tauStarTest with discrete data.", {
  set.seed(238)
  sims = 10
  n = 2000
  for (i in 1:sims) {
    p = runif(sample(1:6, 1))
    p = p / sum(p)
    q = runif(sample(1:6, 1))
    q = q / sum(q)
    x = sample(length(p), n, replace = T, prob = p)
    y = sample(length(q), n, replace = T, prob = q)
    tval = tStar(x, y)
    a = tauStarTest(x, y, mode = "discrete")
    expect_equal(a$x, x)
    expect_equal(a$y, y)
    expect_equal(a$tStar, tval)
    expect_true(abs(a$pVal - (1 - pDisHoeffInd(n * tval, p, q))) <= 10^-2)
    expect_equal(a$mode, "discrete")
  }

  p = runif(sample(1:6, 1))
  p = p / sum(p)
  q = runif(sample(1:6, 1))
  q = q / sum(q)
  x = sample(length(p), n, replace = T, prob = p)
  y = sample(length(q), n, replace = T, prob = q)
  tval = tStar(x, y)
  a = tauStarTest(x, y)
  expect_equal(a$x, x)
  expect_equal(a$y, y)
  expect_equal(a$tStar, tval)
  expect_true(abs(a$pVal - (1 - pDisHoeffInd(n * tval, p, q))) <= 10^-2)
  expect_equal(a$mode, "auto-discrete")

  expect_error(tauStarTest(rnorm(10), rnorm(9), mode="discrete"))
})

test_that("tauStarTest with mixed data.", {
  set.seed(238)
  sims = 10
  n = 200
  for (i in 1:sims) {
    p = runif(sample(1:6, 1))
    p = p / sum(p)
    if (sample(c(T,F), 1)) {
      x = sample(length(p), n, replace = T, prob = p)
      y = rnorm(n)
    } else {
      y = sample(length(p), n, replace = T, prob = p)
      x = rnorm(n)
    }
    tval = tStar(x, y)
    a = tauStarTest(x, y, mode = "mixed")
    expect_equal(a$x, x)
    expect_equal(a$y, y)
    expect_equal(a$tStar, tval)
    expect_true(abs(a$pVal - (1 - pMixHoeffInd(n * tval, p))) <= 10^-2)
    expect_equal(a$mode, "mixed")
  }

  p = runif(sample(1:6, 1))
  p = p / sum(p)
  x = sample(length(p), n, replace = T, prob = p)
  y = rnorm(n)
  tval = tStar(x, y)
  a = tauStarTest(x, y, mode = "auto")
  expect_equal(a$x, x)
  expect_equal(a$y, y)
  expect_equal(a$tStar, tval)
  expect_true(abs(a$pVal - (1 - pMixHoeffInd(n * tval, p))) <= 10^-2)
  expect_equal(a$mode, "auto-mixed")

  x = rnorm(n)
  expect_warning(tauStarTest(x, y, mode = "mixed"))

  expect_error(tauStarTest(rnorm(10), rnorm(9), mode = "mixed"))
  expect_error(tauStarTest(c(1,2,3,4,1), c(1,2,3,4,1), mode = "mixed"))
})

test_that("tauStarTest with a permutation test.", {
  set.seed(238)
  sims = 10
  n = 50
  resamples = 123
  for (i in 1:sims) {
    x = rnorm(n)
    y = rnorm(n)
    tval = tStar(x, y)
    a = tauStarTest(x, y, mode = "permutation", resamples = resamples)
    expect_equal(a$x, x)
    expect_equal(a$y, y)
    expect_equal(a$tStar, tval)
    expect_true(!is.null(a$pVal))
    expect_equal(a$resamples, resamples)
    expect_equal(a$mode, "permutation")
  }

  expect_error(tauStarTest(rnorm(10), rnorm(9), mode = "permutation"))
})
