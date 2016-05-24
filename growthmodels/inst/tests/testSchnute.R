##
##  Tests for the Schnute growth model
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

context("Schnute growth model")

MAXERROR <- 1e-6

test_that("Schnute growth model values", {
  expected   <- c(1.002476, 1.011048,  1.048606,  1.202606,   1.732051,
                  3.156482, 6.416469, 13.454897, 28.422836, 570.535348)
  parameters <- c(1, 2, 3, 0.5)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 4.0)
  
  expect_that(schnute(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(1.643888, 1.644360, 1.646471, 1.655800,  1.695218,
                  1.837775, 2.205393, 2.862040, 3.825090, 12.662398)
  parameters <- c(12, 2, 3, .20)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 4.0)
  
  expect_that(schnute(time, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Schnute growth model values", {
  parameters <- c(1, 2, 3, 0.5)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 4.0)
  size       <- schnute(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(schnute.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 2, 3, .20)
  size       <- schnute(time, parameters[1], parameters[2], parameters[3], parameters[4])
  
  expect_that(schnute.inverse(size, parameters[1], parameters[2], parameters[3], parameters[4]),
              equals(time, tolerance = MAXERROR))
})
