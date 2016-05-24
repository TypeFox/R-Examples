library(testthat)
library(ff)

context("sum")

test_that("Sum ff works",{
  x <- runif(100)
  fx <- ff(x)
  expect_equal(sum(x), sum(fx))
})

test_that("Sum ff works with na=TRUE works",{
  x <- 1
  fx <- ff(x)
  expect_equal(sum(x, na.rm=TRUE), sum(fx, na.rm=TRUE))
})