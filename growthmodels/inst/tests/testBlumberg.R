##
##  Tests for the Blumberg growth model
##
##  Created by Daniel Rodríguez Pérez on 14/9/2013.
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

context("Blumberg growth model")

MAXERROR <- 1e-6

test_that("Blumberg growth model values", {
  expected   <- c(0.00000000, 0.05882353, 0.33333333, 0.62790698, 0.80000000,
                  0.88652482, 0.93103448, 0.95543175, 0.96969697, 0.97852349)
  parameters <- c(1, 2, 3)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(blumberg(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(0.000000, 1.292630, 1.666667, 1.923072, 2.122552,
                  2.287412, 2.428652, 2.552605, 2.663284, 2.763410)
  parameters <- c(10, 5, .43)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(blumberg(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(0.8000000, 0.8865248, 0.9310345, 0.9554318, 0.9696970,
                  0.9785235, 0.9842520, 0.9881218, 0.9908257, 0.9927700)
  parameters <- c(1, 2, 3, 2)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(blumberg(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(2.122552, 2.287412, 2.428652, 2.552605, 2.663284,
                  2.763410, 2.854921, 2.939251, 3.017494, 3.090503)
  parameters <- c(10, 5, .43, 2)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(blumberg(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Blumberg growth model values", {
  parameters <- c(1, 2, 3)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  size       <- blumberg(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(blumberg.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(10, 5, .43)
  size       <- blumberg(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(blumberg.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(1, 2, 3, 2)
  size       <- blumberg(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(blumberg.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(10, 5, .43, 2)
  size       <- blumberg(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(blumberg.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})
