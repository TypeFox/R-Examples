context("imap iterator")

pow <- function(x, y) {
  x^y
}

test_that("imap over two numeric vectors of equal length", {
  it <- imap(pow, c(2, 3, 10), c(5, 2, 3))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})

test_that("imap over two numeric vectors of unequal length", {
  it <- imap(pow, c(2, 3, 10), c(5, 2, 3, 42))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})

test_that("imap over two lists of equal length", {
  it <- imap(pow, list(2, 3, 10), list(5, 2, 3))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})

test_that("imap over two lists of unequal length", {
  it <- imap(pow, list(2, 3, 10), list(5, 2, 3, 42))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})
