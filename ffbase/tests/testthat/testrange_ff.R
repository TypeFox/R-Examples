library(testthat)
library(ff)

context("range")

test_that("Range ff works",{
  x <- runif(100)
  fx <- ff(x)
  expect_equal(range(x), range(fx))
  
  is.na(x) <- sample(100, 10)
  fx <- ff(x)
  expect_equal(range(x), range(fx))
  expect_equal(range(x, na.rm=TRUE), range(fx, na.rm=TRUE))
})