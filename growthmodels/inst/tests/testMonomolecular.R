##
##  Tests for the monomolecular growth model
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

context("Monomolecular growth model")

MAXERROR <- 1e-6

test_that("Monomolecular growth model values", {
  expected   <- c(-8.058576e+002, -1.790343e+002, -3.917107e+001,
                  -7.963378e+000, -1, 5.537397e-001, 9.004259e-001,
                  9.777820e-001, 9.950425e-001)
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(monomolecular(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
  
  expected   <- c(-6.431778e+002, -2.290264e+002, -7.666867e+001,
                  -2.061938e+001, 0, 7.585447e+000, 1.037598e+001,
                  1.140256e+001, 1.178021e+001)
  parameters <- c(12, 1, 2)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  
  expect_that(monomolecular(time, parameters[1], parameters[2], parameters[3]),
              equals(expected, tolerance = MAXERROR))
})

test_that("Inverse Monomolecular growth model values", {
  parameters <- c(1, 2, 3)
  time       <- c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
  size       <- monomolecular(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(monomolecular.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
  
  parameters <- c(12, 1, 2)
  size       <- monomolecular(time, parameters[1], parameters[2], parameters[3])
  
  expect_that(monomolecular.inverse(size, parameters[1], parameters[2], parameters[3]),
              equals(time, tolerance = MAXERROR))
})
