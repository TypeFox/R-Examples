##
##  Tests for the Log-logistic growth model
##
##  Created by Daniel Rodríguez Pérez on 30/8/2013.
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

context("Log-logistic growth model")

MAXERROR <- 1e-6

test_that("Log-logistic growth model invalid values", {
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5)
  
  result <- loglogistic(time, parameters[1], parameters[2], parameters[3])
  
  for (i in 1:length(time)) {
    expect_that(is.na(result[i]), is_true())
  }
})

test_that("Log-logistic growth model in t = 0", {
  expect_that(loglogistic(0,  1, 2, 3),   equals(0, tolerance = MAXERROR))
  expect_that(loglogistic(0, 10, 1, 0.3), equals(0, tolerance = MAXERROR))
  expect_that(loglogistic(0, 12, 2, 0.5), equals(0, tolerance = MAXERROR))
})

test_that("Log-logistic growth model values", {
  expected   <- c(0.00000000, 0.05882353, 0.33333333, 0.62790698, 0.80000000,
                  0.93103448, 0.96969697, 0.98425197, 0.99082569, 0.99420290)
  parameters <- c(1, 2, 3)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
  
  expect_that(loglogistic(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c( 0.000000,  2.400000,  6.000000,  8.307692,  9.600000,
                  10.800000, 11.294118, 11.538462, 11.675676, 11.760000)
  parameters <- c(12, 1, 2)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
  
  expect_that(loglogistic(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Log-logistic growth model values", {
  parameters <- c(1, 2, 3)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
  size       <- loglogistic(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(loglogistic.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 1, 2)
  size       <- loglogistic(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(loglogistic.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
})