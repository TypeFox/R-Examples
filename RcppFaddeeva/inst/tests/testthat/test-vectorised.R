library(RcppFaddeeva)
library(testthat)

z <- 1:10 + 1i

context("Checking the vectorised interface")

test_that("All functions are vectorised", {
  
  expect_equal(length(erfcx(z)), 10)
  expect_equal(length(erf(z)), 10)
  expect_equal(length(erfi(z)), 10)
  expect_equal(length(erfc(z)), 10)
  expect_equal(length(Dawson(z)), 10)
  expect_equal(length(Faddeeva_w(z)), 10)
                      
})
