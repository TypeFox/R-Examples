#unit testing
library(testthat)
library(mbbefd)

context("Check some values with Mahler")

test_that("some values", {
  expect_equal(round(ecmbbefd(x=1/2,a = .2,b = .04),2), 0.68)
  expect_equal(round(dmbbefd(x = .1, a=.2, b=.04),2), .65)
  expect_equal(round(1-pmbbefd(q= .6, a=.2, b=.04),4),.5043)  
})


test_that("mbbefdExposure", {
  expect_equal(ecmbbefd(0.5, a=1, b=1), 0.5)  
  expect_equal(ecmbbefd(0.5, a=0, b=1), 0.5) 
  expect_equal(ecmbbefd(0.5, a=0, b=0), NaN)  
  expect_equal(ecmbbefd(0.5, a=1, b=0), NaN)  
  
  expect_equal(ecMBBEFD(0.5, g=1, b=1), 0.5)  
  expect_equal(ecMBBEFD(0.5, g=1, b=0), 0.5)
  
})

