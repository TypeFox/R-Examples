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

context("nsprcomp.pca")

test_that("PCA equivalence, square", {
    set.seed(1)
    d <- 10
    X <- matrix(rnorm(d*d), d)
    nspc.model <- nsprcomp(X, ncomp = d-1, em_tol = 1e-10)
    pc.model <- prcomp(X)
    
    rot_nrm <- norm(abs(nspc.model$rotation) - abs(pc.model$rotation[ ,1:(d-1)]), "F")
    expect_true(rot_nrm < 1e-3)
    x_nrm <- norm(abs(nspc.model$x) - abs(pc.model$x[ ,1:(d-1)]), "F")
    expect_true(x_nrm < 1e-3)
    sdev_nrm <- sqrt(sum((nspc.model$sdev - pc.model$sdev[1:(d-1)])^2))
    expect_true(sdev_nrm < 1e-3)
})

test_that("PCA equivalence, fat", {
    set.seed(1)
    d <- 10
    n <- 5
    X <- matrix(rnorm(n*d), n)
    nspc.model <- nsprcomp(X, ncomp = n-1, em_tol = 1e-10)
    pc.model <- prcomp(X)
    
    rot_nrm <- norm(abs(nspc.model$rotation) - abs(pc.model$rotation[ ,1:(n-1)]), "F")
    expect_true(rot_nrm < 1e-3)
    x_nrm <- norm(abs(nspc.model$x) - abs(pc.model$x[ ,1:(n-1)]), "F")
    expect_true(x_nrm < 1e-3)
    sdev_nrm <- sqrt(sum((nspc.model$sdev - pc.model$sdev[1:(n-1)])^2))
    expect_true(sdev_nrm < 1e-3)
})

test_that("PCA equivalence, skinny", {
    set.seed(1)
    d <- 5
    n <- 10
    X <- matrix(rnorm(n*d), n)
    nspc.model <- nsprcomp(X, ncomp = d-1, em_tol = 1e-10)
    pc.model <- prcomp(X)
    
    rot_nrm <- norm(abs(nspc.model$rotation) - abs(pc.model$rotation[ ,1:(d-1)]), "F")
    expect_true(rot_nrm < 1e-3)
    x_nrm <- norm(abs(nspc.model$x) - abs(pc.model$x[ ,1:(d-1)]), "F")
    expect_true(x_nrm < 1e-3)
    sdev_nrm <- sqrt(sum((nspc.model$sdev - pc.model$sdev[1:(d-1)])^2))
    expect_true(sdev_nrm < 1e-3)
})


test_that("sdev tolerance early stopping", {
    set.seed(1)
    X <- matrix(runif(10*10), 10)
    
    nspc <- nsprcomp(X, tol = 0.3)
    ncomp <- length(nspc$sdev)
    expect_true(nspc$sdev[ncomp]/nspc$sdev[1] >= 0.3)
    expect_true(ncol(nspc$rotation) == ncomp)
    expect_true(ncol(nspc$q) == ncomp)
})

test_that("rank of matrix smaller than ncomp", {
    a <- 1:5
    X <- a %o% a
    
    nspc <- nsprcomp(X, ncomp = 3)
    expect_true(length(nspc$sdev) == 1)
    expect_true(ncol(nspc$rotation) == 1)
    expect_true(ncol(nspc$q) == 1)
})

test_that("integer weighted PCA equal to repeated observations", {
    set.seed(1)
    X <- matrix(runif(9), 3)
    Y <- rbind(X, X[3,])
    
    set.seed(1)
    nspc <- nsprcomp(Y, center = FALSE)
    nspc.weighted <- nsprcomp(X, omega = c(1,1,2), center = FALSE)
    w1 <- nspc$rotation[ ,1]
    w2 <- nspc.weighted$rotation[ ,1]
    expect_true(sum(abs(w1-w2)) < 1e-3)
})

test_that("reconstruction", {
    set.seed(1)
    X <- matrix(runif(5*5), 5)
    nspc <- nsprcomp(X)
    X_hat <- predict(nspc)%*%ginv(nspc$rotation) + matrix(1,5,1) %*% nspc$center

    expect_true(norm(X - X_hat, type="F") < 1e-3)
})

test_that("weighted approximation error", {
    set.seed(1)
    X <- scale(matrix(runif(5*5), 5))
    nspc <- nsprcomp(X, omega = c(1,1,1,1,5), ncomp = 2)
    X_hat <- predict(nspc)%*%ginv(nspc$rotation)
    
    nrm <- rowSums((X - X_hat)^2)
    expect_true(which.min(nrm) == 5)
})

test_that("sequential component computation", {
    set.seed(1)
    d <- 5
    X <- scale(matrix(runif(d*d), d))
    
    pc.model <- prcomp(X)
    nspc.model <- NULL
    for (pp in seq(d-1)) {
        nspc.model <- nsprcomp(X, ncomp = pp, em_tol = 1e-10, 
                               partial_model = nspc.model)
    }
    
    rot_nrm <- norm(abs(nspc.model$rotation) - abs(pc.model$rotation[ ,1:(d-1)]), "F")
    expect_true(rot_nrm < 1e-3)
    x_nrm <- norm(abs(nspc.model$x) - abs(pc.model$x[ ,1:(d-1)]), "F")
    expect_true(x_nrm < 1e-3)
    sdev_nrm <- sqrt(sum((nspc.model$sdev - pc.model$sdev[1:(d-1)])^2))
    expect_true(sdev_nrm < 1e-3)
})

test_that("continuing an early stopped model", {
    set.seed(1)
    d <- 10
    X <- scale(matrix(runif(d*d), d))
    
    pc.model <- prcomp(X)
    nspc.model <- nsprcomp(X, tol = 0.5, em_tol = 1e-10)
    nspc.model <- nsprcomp(X, em_tol = 1e-10, partial_model = nspc.model)
    
    rot_nrm <- norm(abs(nspc.model$rotation) - abs(pc.model$rotation[ ,1:(d-1)]), "F")
    expect_true(rot_nrm < 1e-3)
    x_nrm <- norm(abs(nspc.model$x) - abs(pc.model$x[ ,1:(d-1)]), "F")
    expect_true(x_nrm < 1e-3)
    sdev_nrm <- sqrt(sum((nspc.model$sdev - pc.model$sdev[1:(d-1)])^2))
    expect_true(sdev_nrm < 1e-3)
})
