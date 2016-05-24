context("itabulate iterator")

test_that("itabulate functions properly with default values", {
  it <- itabulate(f=function(x) x + 1)
  expect_equal(take(it, 4), as.list(2:5))
})

test_that("itabulate functions properly with a given start value", {
  it <- itabulate(f=function(x) x^2, start=-3)
  expect_equal(take(it, 6), list(9, 4, 1, 0, 1, 4))
})

test_that("itabulate functions properly with given start and step values", {
  it <- itabulate(abs, start=-5, step=2)
  expect_equal(take(it, 6), list(5, 3, 1, 1, 3, 5))
})

test_that("itabulate functions properly with a decreasing sequence", {
  it <- itabulate(exp, start=6, step=-2)
  expect_equal(take(it, 4), as.list(exp(seq(6, 0, by=-2))))
})
