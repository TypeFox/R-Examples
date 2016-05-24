##
##  Tests for the Chapman-Richards growth model
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

context("Chapman-Richards growth model")

MAXERROR <- 1e-6

test_that("Chapman-Richards growth model invalid values", {
  parameters <- c(1, 2, 3, 4)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0)
  
  result <- chapmanRichards(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  for (i in 1:length(time)) {
    expect_that(is.na(result[i]), is_true())
  }
})

test_that("Chapman-Richards growth model values", {
  expected   <- c(1.217769e+000, 1.035581e+000, 1.007518e+000, 1.001658e+000,
                  1.000369e+000, 1.000082e+000, 1.000018e+000, 1.000004e+000,
                  1.000001e+000)
  parameters <- c(1, 2, 3, 4)
  time       <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  
  expect_that(chapmanRichards(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Chapman-Richards growth model values", {
  parameters <- c(1, 2, 3, 4)
  time       <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5)
  size       <- chapmanRichards(time, parameters[1], parameters[2],
                                parameters[3], parameters[4])
  
  expect_that(chapmanRichards.inverse(size, parameters[1], parameters[2],
                                      parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(1, 1, 2, 2)
  size       <- chapmanRichards(time, parameters[1], parameters[2],
                                parameters[3], parameters[4])
  
  expect_that(chapmanRichards.inverse(size, parameters[1], parameters[2],
                                      parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})
