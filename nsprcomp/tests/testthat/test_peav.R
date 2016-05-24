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

context("percentage explained additional variance")

test_that("PCA sums to 1", {
    set.seed(1)
    
    X <- matrix(runif(10*5), 10)
    pc <- prcomp(X)
    expect_equal(sum(peav(X, pc$rotation)), 1)
    
    X <- matrix(runif(10*5), 5)
    pc <- prcomp(X)
    expect_equal(sum(peav(X, pc$rotation)), 1)
})

test_that("arbitrary rotation sums to 1", {
    set.seed(1)
    
    X <- matrix(runif(10*5), 10)
    W <- qr.Q(qr(matrix(rnorm(5*5),5)))
    expect_equal(sum(peav(X, W)), 1)
    
    X <- matrix(runif(10*5), 5)
    W <- qr.Q(qr(matrix(rnorm(10*10),10)))
    expect_equal(sum(peav(X, W)), 1)
})

