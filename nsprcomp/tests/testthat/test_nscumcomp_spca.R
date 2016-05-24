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

context("nscumcomp.spca")

test_that("cardinality", {
    set.seed(1)
    X <- matrix(rnorm(20*10), 20)
    
    nscc <- nscumcomp(X, ncomp = 5, gamma = 1)
    expect_equal(sum(cardinality(nscc$rotation)), 50)
    
    nscc <- nscumcomp(X, ncomp = 1, gamma = 1, k = 5)
    expect_equal(sum(cardinality(nscc$rotation)), 5)
    
    nscc <- nscumcomp(X, ncomp = 5, gamma = 50, k = 10)
    expect_equal(sum(cardinality(nscc$rotation)), 10)
})

test_that("reconstruction", {
    set.seed(1)
    X <- matrix(runif(5*5), 5)
    nscc <- nscumcomp(X, ncomp = 5, k = 20, gamma = 1)
    X_hat <- predict(nscc)%*%ginv(nscc$rotation) + matrix(1,5,1) %*% nscc$center
    
    expect_true(norm(X - X_hat, type="F") < 1e-3)
})

test_that("weighted approximation error", {
    set.seed(1)
    X <- scale(matrix(runif(5*5), 5))
    nscc <- nscumcomp(X, omega = c(1,1,1,1,5), ncomp = 3, k = 15, gamma = 1)
    X_hat <- predict(nscc)%*%ginv(nscc$rotation)
    
    nrm <- rowSums((X - X_hat)^2)
    expect_true(which.min(nrm) == 5)
})
