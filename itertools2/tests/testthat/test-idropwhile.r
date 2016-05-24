context("idropwhile iterator")

test_that("Apply idropwhile to an integer sequence", {
  not_too_large <- function(x) {
    x <= 3
  }
  it <- idropwhile(not_too_large, 1:10)

  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 5)
  expect_equal(nextElem(it), 6)
  expect_equal(nextElem(it), 7)
  expect_equal(nextElem(it), 8)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 10)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply idropwhile to an integer sequence using anonymous function", {
  it <- idropwhile(function(x) x <= 10, seq(2, 20, by=2))

  expect_equal(nextElem(it), 12)
  expect_equal(nextElem(it), 14)
  expect_equal(nextElem(it), 16)
  expect_equal(nextElem(it), 18)
  expect_equal(nextElem(it), 20)
  expect_error(nextElem(it), "StopIteration")
})
