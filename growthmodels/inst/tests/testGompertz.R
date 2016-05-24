##
##  Tests for the Gompertz growth model
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

context("Gompertz growth model")

MAXERROR <- 1e-6

test_that("Gompertz growth model values", {
  expected   <- c(0, 6.488035e-079, 3.580340e-018, 1.280131e-004, 1.353353e-001,
                  6.400171e-001, 9.052228e-001, 9.780270e-001, 9.950548e-001)
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(gompertz(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(2.330805e-023, 2.270614e-008, 7.415748e-003, 7.918564e-001,
                  4.414553e+000, 8.306408e+000, 1.048108e+001, 1.141718e+001,
                  1.178221e+001)
  parameters <- c(12, 1, 2)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(gompertz(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Gompertz growth model values", {
  parameters <- c(1, 2, 3)
  time       <- c(-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- gompertz(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(gompertz.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 1, 2)
  size       <- gompertz(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(gompertz.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
})
