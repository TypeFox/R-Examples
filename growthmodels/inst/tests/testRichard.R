##
##  Tests for the Richard growth model
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

context("Richard growth model")

MAXERROR <- 1e-6

test_that("Richard growth model invalid values", {
  parameters <- c(1, -2, 3, 4)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0)
  
  result <- richard(time, parameters[1], parameters[2], parameters[3],
                    parameters[4])
  
  for (i in 1:length(time)) {
    expect_that(is.na(result[i]), is_true())
  }
})

test_that("Richard growth model values", {
  expected   <- c(1.875713e-001, 2.726213e-001, 3.947771e-001, 5.628574e-001,
                  7.598357e-001, 9.118815e-001, 9.765486e-001, 9.945214e-001,
                  9.987644e-001)
  parameters <- c(1, 2, 3, 4)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(richard(time, parameters[1], parameters[2], parameters[3],
                      parameters[4]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(1.485410e-002, 6.628579e-002, 2.914668e-001, 1.204411e+000,
                  4, 8.297261e+000, 1.091332e+001, 1.173918e+001, 1.194080e+001)
  parameters <- c(12, 2, 3, 1)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(richard(time, parameters[1], parameters[2], parameters[3],
                      parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Richards growth model values", {
  parameters <- c(1, 2, 3, 4)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- richard(time, parameters[1], parameters[2],
                                parameters[3], parameters[4])
  
  expect_that(richard.inverse(size, parameters[1], parameters[2],
                                      parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 2, 3, 1)
  size       <- richard(time, parameters[1], parameters[2],
                        parameters[3], parameters[4])
  
  expect_that(richard.inverse(size, parameters[1], parameters[2],
                              parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})

test_that("Generalised Richard growth model values", {
  expected   <- c(1.138049e+000, 1.257893e+000, 1.481437e+000, 1.891192e+000,
                  2.519671e+000, 2.925422e+000, 2.993318e+000, 2.999447e+000,
                  2.999955e+000)
  parameters <- c(1, 3, 5, 4, 2, 0)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(generalisedRichard(time, parameters[1], parameters[2], parameters[3],
                                 parameters[4], parameters[5], parameters[6]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(1.200000e+001, 1.199999e+001, 1.199982e+001, 1.199779e+001,
                  1.197314e+001, 1.168460e+001, 9.333333e+000, 5.128151e+000,
                  4.106374e+000)
  parameters <- c(12, 4, 5, 1, 2, 1)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(generalisedRichard(time, parameters[1], parameters[2], parameters[3],
                                 parameters[4], parameters[5], parameters[6]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Richards growth model values", {
  parameters <- c(1, 3, 5, 4, 2, 0)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- generalisedRichard(time, parameters[1], parameters[2],
                                   parameters[3], parameters[4],
                                   parameters[5], parameters[6])
  
  expect_that(generalisedRichard.inverse(size, parameters[1], parameters[2],
                                         parameters[3], parameters[4], parameters[5],
                                         parameters[6]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 4, 5, 1, 2, 1)
  size       <- generalisedRichard(time, parameters[1], parameters[2],
                                   parameters[3], parameters[4],
                                   parameters[5], parameters[6])
  
  expect_that(generalisedRichard.inverse(size, parameters[1], parameters[2],
                              parameters[3], parameters[4], parameters[5],
                                         parameters[6]),
              equals(time, tolerance = MAXERROR))
})
