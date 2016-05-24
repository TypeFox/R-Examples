{###############################################################################
# test-use_method.R
# This file is part of the R package harvestr.
# 
# Copyright 2012 Andrew Redd
# Date: 6/2/2012
# 
# DESCRIPTION
# ===========
# testing for use_method
# 
# LICENSE
# ========
# harvestr is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# dostats. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

library(harvestr)
library(testthat)
context("use_method")

test_that("use_method", {
  exf <- system.file("examples", "use_method.R", package="harvestr")
  if(exf=="") 
    exf <- system.file("inst", "examples", "use_method.R", package="harvestr")
  if(exf!=""){      
      source(exf)
      x <- mr$new(x = 2L)
      x$x <- 2L
      name <- "Hello World!"
      x$name <- "Hello World!"
      
      expect_that(use_method("hello"), is_a("function"))
      expect_that(use_method(hello), is_a("function"))
      expect_that(use_method(hello, 1), is_a("function"))
      expect_that(use_method(hello)(x), is_a("character"))
      expect_that(use_method(times, 3)(x), equals(6))
  } else {expect_that(FALSE, is_true(), "Could not find the example for use_method.")}
})

