##
##  Tests for the von Bertalanffy growth model
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


context("von Bertalanffy growth model")

MAXERROR <- 1e-6

test_that("von Bertalanffy growth model invalid values", {
  parameters <- c(1, 2, 3, 4)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0)
  
  result <- vonBertalanffy(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  for (i in 1:length(time)) {
    expect_that(is.na(result[i]), is_true())
  }
})

test_that("von Bertalanffy growth model values", {
  expected   <- c(1.217769e+000, 1.035581e+000, 1.007518e+000, 1.001658e+000,
                  1.000369e+000, 1.000082e+000, 1.000018e+000, 1.000004e+000,
                  1.000001e+000)
  parameters <- c(1, 2, 3, 4)
  time       <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(vonBertalanffy(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(1.199897e+001, 1.199977e+001, 1.199995e+001, 1.199999e+001,
                  1.200000e+001, 1.200000e+001, 1.200000e+001, 1.200000e+001,
                  1.200000e+001)
  parameters <- c(12, 2, 3, -2)
  time       <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(vonBertalanffy(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Bertalanffy growth model values", {
  parameters <- c(1, 2, 3, 4)
  time       <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  size       <- vonBertalanffy(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(vonBertalanffy.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 2, 3, -2)
  size       <- vonBertalanffy(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(vonBertalanffy.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})
