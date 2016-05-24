context("izip iterator")

test_that("izip properly iterates through two numeric vectors", {
  it <- izip(1:3, 4:6)
  expect_equal(nextElem(it), list(1, 4))
  expect_equal(nextElem(it), list(2, 5))
  expect_equal(nextElem(it), list(3, 6))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip properly iterates through three numeric vectors", {
  it <- izip(1:3, 4:6, 7:9)
  expect_equal(nextElem(it), list(1, 4, 7))
  expect_equal(nextElem(it), list(2, 5, 8))
  expect_equal(nextElem(it), list(3, 6, 9))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip properly iterates through three numeric vectors and one has extra values", {
  it <- izip(1:3, 4:6, 7:10)
  expect_equal(nextElem(it), list(1, 4, 7))
  expect_equal(nextElem(it), list(2, 5, 8))
  expect_equal(nextElem(it), list(3, 6, 9))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip properly iterates through two numeric vectors and a character vector", {
  it <- izip(1:3, 4:10, levels(iris$Species))
  expect_equal(nextElem(it), list(1, 4, "setosa"))
  expect_equal(nextElem(it), list(2, 5, "versicolor"))
  expect_equal(nextElem(it), list(3, 6, "virginica"))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip properly iterates through a numeric vector, a character vector, and a data.frame's columns", {
  it <- izip(1:3, levels(iris$Species), iris)
  expect_equal(nextElem(it), list(1, "setosa", iris[, 1]))
  expect_equal(nextElem(it), list(2, "versicolor", iris[, 2]))
  expect_equal(nextElem(it), list(3, "virginica", iris[, 3]))
  expect_error(nextElem(it), "StopIteration")
})
