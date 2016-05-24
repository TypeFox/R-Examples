##
##  Tests for the Logistic growth model
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

context("Logistic growth model")

MAXERROR <- 1e-6

test_that("Logistic growth model values", {
  expected   <- c(1.237842e-003, 5.523816e-003, 2.428890e-002, 1.003676e-001,
                  3.333333e-001, 6.914385e-001, 9.094430e-001, 9.782649e-001,
                  9.950670e-001 )
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(logistic(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(2.158345e-001, 5.691105e-001, 1.430435e+000, 3.227297e+000,
                  6, 8.772703e+000, 1.056956e+001, 1.143089e+001, 1.178417e+001)
  parameters <- c(12, 1, 2)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(logistic(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Logistic growth model values", {
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- logistic(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(logistic.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 1, 2)
  size       <- logistic(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(logistic.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
})
  
test_that("Generalised Logistic growth model values", {
  expected   <- c(1.000045e+000, 1.000553e+000, 1.006715e+000, 1.078849e+000,
                  1.666667e+000, 2.717962e+000, 2.973407e+000, 2.997790e+000,
                  2.999818e+000)
  parameters <- c(1, 3, 5, 2, 0)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(generalisedLogistic(time, parameters[1], parameters[2], parameters[3], parameters[4], parameters[5]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(1.200000e+001, 1.199998e+001, 1.199980e+001, 1.199751e+001,
                  1.196978e+001, 1.164518e+001, 9, 4.269170e+000, 3.119670e+000)
  parameters <- c(12, 3, 5, 2, 1)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(generalisedLogistic(time, parameters[1], parameters[2], parameters[3], parameters[4], parameters[5]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Generalised Logistic growth model values", {
  parameters <- c(1, 3, 5, 2, 0)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- generalisedLogistic(time, parameters[1], parameters[2], parameters[3],
                                    parameters[4], parameters[5])
  
  expect_that(generalisedLogistic.inverse(size, parameters[1], parameters[2],
                                          parameters[3], parameters[4], parameters[5]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 3, 5, 2, 1)
  size       <- generalisedLogistic(time, parameters[1], parameters[2], parameters[3],
                                    parameters[4], parameters[5])
  
  expect_that(generalisedLogistic.inverse(size, parameters[1], parameters[2],
                                          parameters[3], parameters[4], parameters[5]),
              equals(time, tolerance = MAXERROR))
})
