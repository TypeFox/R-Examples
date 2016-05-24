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

context("additional correlation")

test_that("cancor correlation equivalence", {
    set.seed(1)
    
    equiv <- function (n,d) {
      X <- matrix(runif(n*d), n)
      Y <- matrix(runif(n*d), n)
      cc <- cancor(X, Y)
      expect_equal(acor(X, cc$xcoef, Y, cc$ycoef)$cor, cc$cor)  
    }
    equiv(20, 5)
    equiv(10, 9)
})

test_that("sparse CCA equivalence", {
  set.seed(1)
  
  equiv <- function (n,d) {
    X <- matrix(runif(n*d), n)
    Y <- matrix(runif(n*d), n)
    xpredict <- function(Y, x, cc) {
      en <- glmnet(Y, x, alpha=0.5, intercept=FALSE, dfmax=2)
      V <- coef(en)
      return(V[2:nrow(V), ncol(V)])
    }
    scc <- nscancor(X, Y, xpredict=xpredict, ypredict=xpredict)
    sacc <- acor(X, scc$xcoef, Y, scc$ycoef)
    expect_equal(sacc$cor, scc$cor)  
    expect_equal(sacc$xcoef, scc$xcoef)
    expect_equal(sacc$ycoef, scc$ycoef)
    expect_equal(sacc$xcenter, scc$xcenter)
    expect_equal(sacc$xscale, scc$xscale)
    expect_equal(sacc$ycenter, scc$ycenter)
    expect_equal(sacc$yscale, scc$yscale)
    expect_equal(sacc$xp, scc$xp)
    expect_equal(sacc$yp, scc$yp)
  }
  equiv(10, 5)
  equiv(10, 10)
  equiv(5, 10)
})

test_that("non-negative sparse CCA correlation equivalence", {
  set.seed(1)
  
  equiv <- function (n,d) {
    X <- matrix(runif(n*d), n)
    Y <- matrix(runif(n*d), n)
    xpredict <- function(Y, x, cc) {
      en <- glmnet(Y, x, alpha=0.5, intercept=FALSE, dfmax=2, lower.limits=0)
      V <- coef(en)
      return(V[2:nrow(V), ncol(V)])
    }
    nscc <- nscancor(X, Y, xpredict=xpredict, ypredict=xpredict)
    expect_equal(acor(X, nscc$xcoef, Y, nscc$ycoef)$cor, nscc$cor)  
  }
  equiv(10, 5)
  equiv(10, 10)
  equiv(5, 10)
})