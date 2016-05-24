context("itakewhile iterator")

test_that("Apply itakewhile to an integer sequence", {
  not_too_large <- function(x) {
    x <= 5
  }
  it <- itakewhile(not_too_large, 1:100)

  expect_equal(nextElem(it), 1)
  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 3)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 5)
  expect_error(nextElem(it), "StopIteration")
})

test_that("Apply itakewhile to an integer sequence using anonymous function", {
  it <- itakewhile(function(x) x <= 10, seq(2, 100, by=2))

  expect_equal(nextElem(it), 2)
  expect_equal(nextElem(it), 4)
  expect_equal(nextElem(it), 6)
  expect_equal(nextElem(it), 8)
  expect_equal(nextElem(it), 10)
  expect_error(nextElem(it), "StopIteration")
})
