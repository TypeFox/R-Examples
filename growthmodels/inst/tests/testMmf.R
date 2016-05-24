##
##  Tests for the Morgan-Mercer-Flodin growth model
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

context("Morgan-Mercer-Flodin growth model")

MAXERROR <- 1e-6

test_that("Morgan-Mercer-Flodin growth model values", {
  expected   <- c(1.157895, 1.372093, 1.750000, 1.979592, 2.000000,
                  1.979592, 1.750000, 1.372093, 1.157895, 1.035714)
  parameters <- c(1, 2, 3, 4)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0)
  
  expect_that(mmf(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(-18.000000, -8.000000, -3.000000, 0.000000, 2.000000,
                    3.428571,  4.500000,  5.333333, 6.000000, 7.000000)
  parameters <- c(12, 2, 3, 1)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0)
  
  expect_that(mmf(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Morgan-Mercer-Flodin growth model values", {
  parameters <- c(1, 2, 3, 4)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0)
  size       <- mmf(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(mmf.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 2, 3, 1)
  size       <- mmf(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(mmf.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})
