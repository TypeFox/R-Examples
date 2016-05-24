context("istarmap iterator")

pow <- function(x, y) {
  x^y
}

test_that("istarmap over a list of two numeric vectors of equal length", {
  it <- istarmap(pow, list(c(2, 3, 10), c(5, 2, 3)))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})

test_that("istarmap over a list of two numeric vectors of unequal length", {
  it <- istarmap(pow, list(c(2, 3, 10), c(5, 2, 3, 42)))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})

test_that("istarmap over a list of two lists of equal length", {
  it <- istarmap(pow, list(list(2, 3, 10), list(5, 2, 3)))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})

test_that("istarmap over a list of two lists of unequal length", {
  it <- istarmap(pow, list(list(2, 3, 10), list(5, 2, 3, 42)))

  expect_equal(nextElem(it), 32)
  expect_equal(nextElem(it), 9)
  expect_equal(nextElem(it), 1000)

  expect_error(nextElem(it), "StopIteration")
})

test_that("istarmap over the rows of a data.frame", {
  # Computes sum of each row in the iris data set
  it <- istarmap(sum, iris[, -5])

  expect_equal(nextElem(it), sum(iris[1, -5]))
  expect_equal(nextElem(it), sum(iris[2, -5]))
  expect_equal(nextElem(it), sum(iris[3, -5]))
  expect_equal(nextElem(it), sum(iris[4, -5]))
  expect_equal(nextElem(it), sum(iris[5, -5]))
  expect_equal(nextElem(it), sum(iris[6, -5]))
  # and so on...
})

