
##

context("Utility functions -- decibel conversions")

test_that("conversions are accurate",{
  
  expect_equal(dB(2), 10*log10(2))
  expect_equal(dB(1), 0)
  # from dB to amp or power
  expect_equal(dB(0, invert=TRUE), 1)
  expect_equal(dB(0, invert=TRUE, is.power = TRUE), 1)
  # power ratio
  expect_equal(dB(70, invert = TRUE, is.power = FALSE), 1e7)
  # amplitude ratio
  expect_equal(round(dB(70, invert = TRUE, is.power = TRUE)), 3162)
})

##

context("Utility functions -- variance of difference series")

test_that("difference series are accurate",{
  X <- 1:10
  expect_equal(vardiff(X), 0)
  expect_equal(varddiff(X), 0)
})

##

context("Utility functions -- matrix generation/manipulation")

test_that("colvec reshapes correctly", {
  m <- na_mat(2,2)
  expect_is(colvec(m), 'matrix')
  expect_equal(dim(colvec(m)), c(4,1))
})

test_that("rowvec reshapes correctly", {
  m <- na_mat(2,2)
  expect_is(rowvec(m), 'matrix')
  expect_equal(dim(rowvec(m)), c(1,4))
})

test_that("NA matrices are assembled correctly", {
  
  expect_is(na_mat(1),'matrix')
  expect_error(na_mat(-1))
  expect_error(na_mat(1,-1))
  
  expect_equal(dim(na_mat(2)), c(2,1))
  expect_equal(dim(na_mat(2,2)), c(2,2))
  expect_equal(dim(na.omit(na_mat(2,2))), c(0,2))
  
  expect_equal(range(na_mat(2,2)), c(NA_real_, NA_real_))
  expect_warning(range(na_mat(2,2), na.rm=TRUE))
               
})
  
test_that("zeros and ones are assembled correctly", {
  
  expect_is(zeros(1),'matrix')
  expect_equal(dim(zeros(2)), c(2,1))
  expect_equal(range(zeros(2)), c(0,0))
  
  expect_is(ones(1),'matrix')
  expect_equal(dim(ones(2)), c(2,1))
  expect_equal(range(ones(2)), c(1,1))
  
})

##

context("Utility functions -- modulo division")

test_that("values are accurate", {

  expect_equal(mod(1,2), 1%%2)
  expect_equal(mod(1+1,2+1), (1+1)%%(2+1))
  
  mx <- 3
  X <- -mx:mx
  
  expect_equal(modulo_floor(X), modulo_floor(X, 2))
  expect_equal(modulo_floor(X,-2), modulo_floor(X, 2))
  
  n <- 1; expect_equal(range(modulo_floor(X, n)), range(X))
  n <- mx; expect_equal(range(modulo_floor(X, n)), n*c(-1,1))
  n <- mx+1; expect_equal(range(modulo_floor(X, n)), c(0,0))
  expect_equal(range(modulo_floor(X, NA)), c(0,0))
  
})

test_that("error checking works", {
  
  expect_warning(modulo_floor(1:10, Inf))
  
  expect_error(modulo_floor(1:10,0))
  expect_error(modulo_floor(1:10,NULL))
  expect_error(modulo_floor(1:10,1:2))
  expect_error(modulo_floor(1:10,-2:-1))
  
})

##
