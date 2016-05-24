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

context("nscancor")

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
    nscc <- nscancor(X, Y, xpredict=xpredict, ypredict=xpredict, 
                     iter_tol=1e-10, iter_max=500)
    
    expect_true(normv(cc$cor - nscc$cor) < 1e-3)
    expect_true(norm(abs(cc$xcoef) - abs(nscc$xcoef), "F") < 1e-3)
    expect_true(norm(abs(cc$ycoef) - abs(nscc$ycoef), "F") < 1e-3)
    expect_true(normv(cc$xcenter - nscc$xcenter) < 1e-3)
    expect_true(normv(cc$ycenter - nscc$ycenter) < 1e-3)
})

test_that("corr tolerance early stopping", {
    set.seed(1)
    d <- 5
    n <- 10
    X <- matrix(rnorm(n*d), n)
    Y <- matrix(rnorm(n*d), n)
    
    xpredict = function(Y, x, cc) {
      return(ginv(Y)%*%x)
    } 
    nscc <- nscancor(X, Y, xpredict=xpredict, ypredict=xpredict, cor_tol=0.3)
    ncomp <- length(nscc$cor)
    
    expect_true(nscc$cor[ncomp]/nscc$cor[1] >= 0.3)
    expect_true(ncol(nscc$xcoef) == ncomp)
    expect_true(ncol(nscc$ycoef) == ncomp)
})
 
test_that("rank of matrix smaller than nvar", {
    a <- 1:5
    X <- a %o% a
    
    xpredict = function(Y, x, cc) {
      return(ginv(Y)%*%x)
    } 
    nscc <- nscancor(X, X, xpredict=xpredict, ypredict=xpredict, nvar = 2)
    expect_true(length(nscc$cor) == 1)
    expect_true(ncol(nscc$xcoef) == 1)
    expect_true(ncol(nscc$ycoef) == 1)
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
  nscc <- NULL
  for (pp in seq_along(cc$cor)) {
    nscc <- nscancor(X, Y, xpredict=xpredict, ypredict=xpredict, nvar=pp,
                     iter_tol=1e-10, iter_max=500, partial_model=nscc)
  }
  
  expect_true(normv(cc$cor - nscc$cor) < 1e-3)
  expect_true(norm(abs(cc$xcoef) - abs(nscc$xcoef), "F") < 1e-3)
  expect_true(norm(abs(cc$ycoef) - abs(nscc$ycoef), "F") < 1e-3)
  expect_true(normv(cc$xcenter - nscc$xcenter) < 1e-3)
  expect_true(normv(cc$ycenter - nscc$ycenter) < 1e-3)
})
