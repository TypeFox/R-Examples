{###############################################################################
# test_wargs.R
# Copyright 2012 Andrew Redd
# Date: 6/1/2012
# 
# DESCRIPTION
# ===========
# unit tests for recombine function
# 
# LICENSE
# ========
# dostats is free software: you can redistribute it and/or modify it under the
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
test_that("testing wargs", {
  expect_true(!is.na(
    wargs(mean,na.rm=T)(c(NA, 1:3))
    ))  
  paste1 <- wargs(paste, sep='-')
  expect_equal(
    paste1(1, 2, 3)
    , "1-2-3")
})
test_that("wargs will convey attributes", {
    a <- list(hello='world')
    f <- function(x=1){print(x)}
    attributes(f) <- a
    g <- wargs(f, x=2)
    expect_that(attributes(g), equals(attributes(f)))
})
