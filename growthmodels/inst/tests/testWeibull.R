##
##  Tests for the Weibull growth model
##
##  Created by Daniel Rodríguez Pérez on 28/7/2013.
##
##  Copyright (c) 2013 Daniel Rodríguez Pérez.
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>
##

context("Weibull growth model")

MAXERROR <- 1e-6

test_that("Weibull growth model values", {
  expected   <- c(1, 9.999995e-001, 9.004259e-001, -6.580582e-001, -1,
                  -6.580582e-001, 9.004259e-001, 9.999995e-001, 1)
  parameters <- c(1, 2, 3, 4)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(weibull(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(-7.948576e+002, -1.680343e+002, -2.817107e+001,
                  3.036622e+000, 10, 1.155374e+001, 1.190043e+001,
                  1.197778e+001, 1.199504e+001)
  parameters <- c(12, 2, 3, 1)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(weibull(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Weibull growth model values", {
  parameters <- c(1, 2, 3, 4)
  time       <- c(0.0, 0.5, 1.0, 1.5)
  size       <- weibull(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(weibull.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 2, 3, 1)
  size       <- weibull(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(weibull.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})

