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

context("multi-domain additional correlation")

test_that("cancor correlation equivalence", {
    set.seed(1)
    
    equiv <- function (n,d) {
      X <- matrix(runif(n*d), n)
      Y <- matrix(runif(n*d), n)
      cc <- cancor(X, Y)
      expect_equal(macor(list(X, Y), list(cc$xcoef, cc$ycoef))$cor[1, 2, ], cc$cor)  
    }
    equiv(20, 5)
    equiv(10, 9)
})

test_that("sparse CCA correlation equivalence", {
  set.seed(1)
  
  equiv <- function (n,d) {
    X <- matrix(runif(n*d), n)
    Y <- matrix(runif(n*d), n)
    xpredict <- function(Y, x, cc) {
      en <- glmnet(Y, x, alpha=0.5, intercept=FALSE, dfmax=2)
      V <- coef(en)
      return(V[2:nrow(V), ncol(V)])
    }
    scc <- nscancor(X, Y, xpredict=xpredict, ypredict=xpredict, nvar=2)
    expect_equal(macor(list(X, Y), list(scc$xcoef, scc$ycoef))$cor[1, 2, ], scc$cor)  
  }
  equiv(10, 5)
  equiv(10, 10)
  equiv(5, 10)
})

test_that("sparse mCCA equivalence", {
  set.seed(1)
  
  equiv <- function (n,d) {
    X <- list(matrix(runif(n*d), n), matrix(runif(n*d), n), matrix(runif(n*d), n))
    pred <- function(Y, x, cc) {
      en <- glmnet(Y, x, alpha=0.5, intercept=FALSE, dfmax=2)
      V <- coef(en)
      return(V[2:nrow(V), ncol(V)])
    }
    predict <- list(pred, pred, pred)
    mcc <- mcancor(X, predict=predict, nvar=2)
    macc <- macor(X, mcc$coef)
    expect_equal(macc$cor, mcc$cor)  
    expect_equal(macc$coef, mcc$coef)
    expect_equal(macc$center, mcc$center)
    expect_equal(macc$scale, mcc$scale)
    expect_equal(macc$xp, mcc$xp)
  }
  equiv(10, 5)
  equiv(10, 10)
  equiv(5, 10)
})