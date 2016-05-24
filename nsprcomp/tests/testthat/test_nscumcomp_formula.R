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

context("nscumcomp.formula")

test_that("subset and na_exclude", {
    d <- 3
    n <- 10
    set.seed(1)
    X <- as.data.frame(matrix(runif(n*d), n))
    X[5,1] <- NA
    colnames(X) <- c("One", "Two", "Three")
    
    nscc <- nscumcomp(~ One + Three, data = X, 1:(n-1), 
                      na.exclude, ncomp = 1, gamma = 10)
    expect_equal(unname(nscc$x[5,1]), NA_integer_)
})
