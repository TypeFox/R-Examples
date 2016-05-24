#  Copyright 2013 Christian Sigg
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

context("mcancor")

test_that("cancor equivalence", {
    set.seed(1)
    d <- 5
    n <- 10
    X <- matrix(rnorm(n*d), n)
    Y <- matrix(rnorm(n*d), n)
    
    cc <- cancor(X, Y)
    xpredict = function(Y, x, cc) {
      return(ginv(Y)%*%x)
    } 
    mcc <- mcancor(list(X, Y), predict=list(xpredict, xpredict), nvar = 5,
                     iter_tol=1e-10, iter_max=500)
    
    expect_true(normv(cc$cor - mcc$cor[1, 2, ]) < 1e-3)
    expect_true(norm(abs(cc$xcoef) - abs(mcc$coef[[1]]), "F") < 1e-3)
    expect_true(norm(abs(cc$ycoef) - abs(mcc$coef[[2]]), "F") < 1e-3)
    expect_true(normv(cc$xcenter - mcc$center[[1]]) < 1e-3)
    expect_true(normv(cc$ycenter - mcc$center[[2]]) < 1e-3)
})

test_that("corr tolerance early stopping", {
    set.seed(1)
    d <- 5
    n <- 10
    X <- list(matrix(rnorm(n*d), n), matrix(rnorm(n*d), n), matrix(rnorm(n*d), n))
    
    pred = function(Y, x, cc) {
      return(ginv(Y)%*%x)
    } 
    
    cor_trsh <- function(cor_tol) {
      mcc <- mcancor(X, predict=list(pred, pred, pred), cor_tol=cor_tol)
      m <- length(X)
      nvar <- ncol(mcc$coef[[1]])
      
      sum_corr_np <- sum(mcc$cor[ , , nvar] - diag(3))/2
      sum_corr_1 <- sum(mcc$cor[ , , 1] - diag(3))/2
      expect_true(sum_corr_np/sum_corr_1 >= cor_tol)  
      expect_true(all(sapply(mcc$coef, ncol) == rep(nvar, m)))
    }
    cor_trsh(0)
    cor_trsh(0.6)
    cor_trsh(1)
})
 
test_that("rank of matrix smaller than nvar", {
    a <- 1:5
    X <- a %o% a
    
    pred = function(Y, x, cc) {
      return(ginv(Y)%*%x)
    } 
    mcc <- mcancor(list(X, X), predict=list(pred, pred), nvar = 2)
    
    expect_true(dim(mcc$cor)[3] == 1)
    expect_true(all(sapply(mcc$coef, ncol) == rep(1, 2)))
})

test_that("sequential variable computation", {
  set.seed(1)
  d <- 5
  n <- 10
  X <- matrix(rnorm(n*d), n)
  Y <- matrix(rnorm(n*d), n)
  
  cc <- cancor(X, Y)
  xpredict = function(Y, x, cc) {
    return(ginv(Y)%*%x)
  } 
  mcc <- NULL
  for (pp in seq_along(cc$cor)) {
    mcc <- mcancor(list(X, Y), predict=list(xpredict, xpredict), nvar = pp,
                   iter_tol=1e-10, iter_max=500, partial_model=mcc)
  }
  
  expect_true(normv(cc$cor - mcc$cor[1, 2, ]) < 1e-3)
  expect_true(norm(abs(cc$xcoef) - abs(mcc$coef[[1]]), "F") < 1e-3)
  expect_true(norm(abs(cc$ycoef) - abs(mcc$coef[[2]]), "F") < 1e-3)
  expect_true(normv(cc$xcenter - mcc$center[[1]]) < 1e-3)
  expect_true(normv(cc$ycenter - mcc$center[[2]]) < 1e-3)
})

