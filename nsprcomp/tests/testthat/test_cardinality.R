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

context("cardinality")

test_that("cardinality", {
    A <- rbind(c(1,0),
               c(0,-1),
               c(1,0))
    expect_equal(cardinality(A), c(2,1))
    
    a <- c(1,2,-1,0,1)
    expect_equal(cardinality(a), 4)
    
    b <- rep(0,3)
    expect_equal(cardinality(b), 0)
    
    B <- matrix(rep(0,4), 1)
    expect_equal(cardinality(B), rep(0,4))
})