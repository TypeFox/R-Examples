##
##  Tests for the Brody growth model
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

context("Brody growth model")

MAXERROR <- 1e-6

test_that("Brody growth model values", {
  expected   <- c(2.000000, 1.223130, 1.049787, 1.011109, 1.002479,
                  1.000553, 1.000123, 1.000028, 1.000006, 1.000001)
  parameters <- c(1, 2, 3)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(brody(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(5.000000, 5.967293, 6.747455, 7.376687, 7.884190,
                  8.293511, 8.623646, 8.889914, 9.104669, 9.277879)
  parameters <- c(10, 5, .43)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(brody(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Brody growth model values", {
  parameters <- c(1, 2, 3)
  time       <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  size       <- brody(time, parameters[1], parameters[2], parameters[3])
    
  expect_that(brody.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(10, 5, .43)
  size       <- brody(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(brody.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
})