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

context("additional sdev")

test_that("PCA sdev equivalence", {
    set.seed(1)
    
    X <- matrix(runif(10*5), 10)
    pc <- prcomp(X)
    expect_equal(asdev(X, pc$rotation)$sdev, pc$sdev)
    
    X <- matrix(runif(10*5), 5)
    pc <- prcomp(X)
    expect_equal(asdev(X, pc$rotation)$sdev, pc$sdev)
})

test_that("sparse PCA sdev equivalence", {
    set.seed(1)
    
    X <- matrix(runif(10*5), 10)
    spc <- nsprcomp(X, k=2)
    expect_equal(asdev(X, spc$rotation)$sdev, spc$sdev)
    
    X <- matrix(runif(10*5), 5)
    spc <- nsprcomp(X, k=7)
    expect_equal(asdev(X, spc$rotation)$sdev, spc$sdev)
})

test_that("non-negative sparse PCA sdev equivalence", {
    set.seed(1)
    
    X <- matrix(runif(10*5), 10)
    nspc <- nsprcomp(X, k=2, nneg=TRUE)
    expect_equal(asdev(X, nspc$rotation)$sdev, nspc$sdev)
    
    X <- matrix(runif(10*5), 5)
    nspc <- nsprcomp(X, k=7, nneg=TRUE)
    expect_equal(asdev(X, nspc$rotation)$sdev, nspc$sdev)
})

test_that("total variance of arbitrary rotation", {
    set.seed(1)
    n <- 10
    d <- 5
    
    X <- matrix(rnorm(n*d), n)
    W <- qr.Q(qr(matrix(rnorm(d*d),d)))
    
    expect_equal(sum(asdev(X, W)$sdev^2), sum(scale(X, T, F)^2)/(n-1))
})