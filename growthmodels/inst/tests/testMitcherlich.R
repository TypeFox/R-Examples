##
##  Tests for the Mitcherlich growth model
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

context("Mitcherlich growth model")

MAXERROR <- 1e-6

test_that("Mitcherlich growth model values", {
  expected   <- c(7.777778e-001, 6.150998e-001, 3.333333e-001, -1.547005e-001,
                  -1, -2.464102e+000, -5, -9.392305e+000, -17)
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(mitcherlich(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(1.175000e+001, 1.164645e+001, 1.150000e+001,
                  1.129289e+001, 11, 1.058579e+001, 10, 9.171573e+000, 8)
  parameters <- c(12, 1, 2)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(mitcherlich(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Mitcherlich growth model values", {
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- mitcherlich(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(mitcherlich.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 1, 2)
  size       <- mitcherlich(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(mitcherlich.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
})
