#  Copyright 2012, 2013 Christian Sigg
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

context("nsprcomp.spca")

test_that("cardinality", {
    set.seed(1)
    X <- matrix(rnorm(5*5), 5)
    
    nspc.model <- nsprcomp(X, k = 1)
    expect_true(all(cardinality(nspc.model$rotation) == 1))
    
    nspc.model <- nsprcomp(X, k = 4)
    expect_true(all(cardinality(nspc.model$rotation) <= 4))
    
    nspc.model <- nsprcomp(X, k = 1:5)
    expect_true(all(cardinality(nspc.model$rotation) <= 1:5))
})

test_that("deflation", {
    set.seed(1)
    k <- 4
    d <- 20
    n <- 100
    X = matrix(runif(n*d), n)
    
    nspc <- nsprcomp(X, k = k)
    W <- nspc$rotation
    Xp <- nspc$xp
    for (cc in seq(length(nspc$sdev))) {
        expect_true(sum(abs(Xp%*%W[ ,cc])) < 1e-10)
    }
})

test_that("reconstruction", {
    set.seed(1)
    X <- matrix(runif(5*5), 5)
    nspc <- nsprcomp(X, k = 3)
    X_hat <- predict(nspc)%*%ginv(nspc$rotation) + matrix(1,5,1) %*% nspc$center
    
    expect_true(norm(X - X_hat, type="F") < 1e-3)
})

test_that("weighted approximation error", {
    set.seed(1)
    X <- scale(matrix(runif(5*5), 5))
    nspc <- nsprcomp(X, omega = c(1,1,1,1,5), ncomp = 2, k = 3)
    X_hat <- predict(nspc)%*%ginv(nspc$rotation)
    
    nrm <- rowSums((X - X_hat)^2)
    expect_true(which.min(nrm) == 5)
})