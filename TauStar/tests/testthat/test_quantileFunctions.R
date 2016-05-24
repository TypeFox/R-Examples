library(TauStar)
context("Testing the quantile functions.")

test_that("qHoeffInd function returns correct values", {
  set.seed(9234)
  xVals = rHoeffInd(5)
  for (i in 1:length(xVals)) {
    p = pHoeffInd(xVals[i])
    expect_equal(qHoeffInd(p, 10^-3), xVals[i], tolerance = 10^-3)
  }

  expect_error(qHoeffInd(-.01))
  expect_error(qHoeffInd(1.1))
})

test_that("qDisHoeffInd function returns correct values", {
  set.seed(9234)
  sims = 10
  for (i in 1:sims) {
    p = runif(sample(1:5, 1))
    p = p / sum(p)
    q = runif(sample(1:5, 1))
    q = q / sum(q)
    x = rDisHoeffInd(1, p, q)
    expect_equal(qDisHoeffInd(pDisHoeffInd(x, p, q), p, q, 10^-4),
                 x, tolerance = 10^-3)
  }

  expect_error(qDisHoeffInd(-.01, 1, 1))
  expect_error(qDisHoeffInd(1.1, 1, 1))
})

test_that("qMixHoeffInd function returns correct values", {
  set.seed(9234)
  sims = 10
  for (i in 1:sims) {
    p = runif(sample(1:5, 1))
    p = p / sum(p)
    x = rMixHoeffInd(1, p)
    expect_equal(qMixHoeffInd(pMixHoeffInd(x, p), p, 10^-4),
                 x, tolerance = 10^-3)
  }

  expect_error(qDisHoeffInd(-.01, 1, 1))
  expect_error(qDisHoeffInd(1.1, 1, 1))
})
