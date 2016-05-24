library(testthat)
library(ff)

context("min")

test_that("Min ff works",{
  x <- runif(100) 
  fx <- ff(x)
  expect_equal(min(x), min(fx))
  
  is.na(x) <- sample(100, 10)
  fx <- ff(x)
  expect_equal(min(x), min(fx))
  expect_equal(min(x, na.rm=TRUE), min(fx, na.rm=TRUE))
})

test_that("Min ff works for multiple ff",{
  x <- runif(100)
  y <- runif(20)
  fx <- ff(x)
  fy <- ff(y)
  
  expect_equal(min(x,y), min(fx,fy))
  
  is.na(x) <- sample(100, 10)
  fx <- ff(x)
  expect_equal(min(x), min(fx))
  expect_equal(min(x, na.rm=TRUE), min(fx, na.rm=TRUE))
})

test_that("Min ff works for a range",{
  #x <- runif(100)
  x <- 1:10
  fx <- ff(x)
    
  expect_equal( min(x[2:10])
              , min(fx, range=ri(2,10))
              )  
})
