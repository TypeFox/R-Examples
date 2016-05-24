##
##  Tests for the negative exponential growth model
##
##  Created by Daniel Rodríguez Pérez on 27/7/2013.
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

context("Negative exponential growth model")

MAXERROR <- 1e-6

test_that("Negative exponential growth model values", {
  expected   <- c(-5.359815e+001, -1.908554e+001, -6.389056e+000,
                  -1.718282e+000, 0, 6.321206e-001, 8.646647e-001,
                  9.502129e-001, 9.816844e-001)
  parameters <- c(1, 2)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(negativeExponential(time, parameters[1], parameters[2]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(-7.666867e+001, -4.178027e+001, -2.061938e+001,
                  -7.784655e+000, 0, 4.721632e+000, 7.585447e+000,
                  9.322438e+000, 1.037598e+001)
  parameters <- c(12, 1)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(negativeExponential(time, parameters[1], parameters[2]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse negative exponentia growth model values", {
  parameters <- c(1, 2)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- negativeExponential(time, parameters[1], parameters[2])
  
  expect_that(negativeExponential.inverse(size, parameters[1], parameters[2]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 1)
  size       <- negativeExponential(time, parameters[1], parameters[2])
  
  expect_that(negativeExponential.inverse(size, parameters[1], parameters[2]),
              equals(time, tolerance = MAXERROR))
})
