library(testthat)
library(ff)

context("ffwhich")

test_that("ffwhich on vector works",{
  x <- 10:1
  fx <- ff(x)
  
  idx <- which(x < 5)
  fidx <- ffwhich(fx, fx < 5)
  
  # 
  # dat <- ffdf(x1=x, y1=x)
  # idx <- ffwhich(dat, x1 < 5 & y1 > 2)
  # dat[idx,][,]
  expect_identical(idx, fidx[])   
})

test_that("ffwhich on data.frame works",{
  dat <- data.frame(x1=10:1, y1=10:1)
  fdat <- as.ffdf(dat)
  
   
  fidx <- ffwhich(fdat, x1 < 5 & y1 > 2)
  expect_equivalent(c(7,8), fidx[])   
})

test_that("ffwhich evaluates variables correctly",{
  f <- function(x) {
    a <- 1
    b <- 3
    ffwhich(x, x > a & x < b)
  }
  x <- ff(1:10)
  expect_that(f(x)[], equals(2))
})
