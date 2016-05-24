##
##  Tests for the Stannard growth model
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

context("Stannard growth model")

MAXERROR <- 1e-6

test_that("Stannard growth model values", {
  expected   <- c(0.01831257, 0.08180985, 0.34525776, 0.85501964, 0.99096609,
                  0.99954437, 0.99997730, 0.99999887, 0.99999994, 1.00000000)
  parameters <- c(1, 2, 3, 0.5)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 4.0)
  
  expect_that(stannard(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(0.2192816,   0.9660857,  3.7981941,  9.0791332, 11.5399613,
                  11.9437234, 11.9933643, 11.9992210, 11.9999086, 12.0000000)
  parameters <- c(12, 2, 3, .70)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 4.0)
  
  expect_that(stannard(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Stannard growth model values", {
  parameters <- c(1, 2, 3, 0.5)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 4.0)
  size       <- stannard(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(stannard.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 2, 3, .70)
  size       <- stannard(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(stannard.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})
