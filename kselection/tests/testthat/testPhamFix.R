##
##  Test previous errors on kselection
##
##  Created by Daniel Rodriguez Perez on 10/10/2014.
##
##  Copyright (c) 2014 Daniel Rodriguez Perez.
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

context("Tests previous errors in kselection")

test_that("evaluate data.frames with low rows", {
  x   <- matrix(c(rnorm(5, 2, .1), rnorm(5, 3, .1),
                  rnorm(5, -2, .1), rnorm(5, -3, .1)), 10, 2)
  
  obj <- kselection(x)
  expect_that(class(obj), equals("Kselection"))
  
  expect_warning(kselection(x),
                 "The maximum number of clusters has been reduced from 15 to 9")
})
