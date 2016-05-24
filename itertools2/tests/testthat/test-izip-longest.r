context("izip_longest iterator")

test_that("izip_longest properly iterates through two numeric vectors", {
  it <- izip_longest(1:3, 4:6)
  expect_equal(nextElem(it), list(1, 4))
  expect_equal(nextElem(it), list(2, 5))
  expect_equal(nextElem(it), list(3, 6))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip_longest properly iterates through three numeric vectors", {
  it <- izip_longest(1:3, 4:6, 7:9)
  expect_equal(nextElem(it), list(1, 4, 7))
  expect_equal(nextElem(it), list(2, 5, 8))
  expect_equal(nextElem(it), list(3, 6, 9))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip_longest properly iterates through three numeric vectors and one has extra values", {
  it <- izip_longest(1:3, 4:6, 7:10, fill=NA)
  expect_equal(nextElem(it), list(1, 4, 7))
  expect_equal(nextElem(it), list(2, 5, 8))
  expect_equal(nextElem(it), list(3, 6, 9))
  expect_equal(nextElem(it), list(NA, NA, 10))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip_longest properly iterates through two numeric vectors and a character vector", {
  it <- izip_longest(1:3, 4:7, levels(iris$Species), fill="yooks")
  expect_equal(nextElem(it), list(1, 4, "setosa"))
  expect_equal(nextElem(it), list(2, 5, "versicolor"))
  expect_equal(nextElem(it), list(3, 6, "virginica"))
  expect_equal(nextElem(it), list("yooks", 7, "yooks"))
  expect_error(nextElem(it), "StopIteration")
})

test_that("izip_longest properly iterates through a numeric vector, a character vector, and a data.frame's columns", {
  it <- izip_longest(1:3, levels(iris$Species), iris, fill="zooks")
  expect_equal(nextElem(it), list(1, "setosa", iris[, 1]))
  expect_equal(nextElem(it), list(2, "versicolor", iris[, 2]))
  expect_equal(nextElem(it), list(3, "virginica", iris[, 3]))
  expect_equal(nextElem(it), list("zooks", "zooks", iris[, 4]))
  expect_equal(nextElem(it), list("zooks", "zooks", iris[, 5]))
  expect_error(nextElem(it), "StopIteration")
})
