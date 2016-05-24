{###############################################################################
# test_make_new_id.R
# Copyright 2012 Andrew Redd
# This file is part of the R package dostats
# Date: 5/30/2012
# 
# DESCRIPTION
# ===========
# testing for make_new_id.
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
{## setup
library(testthat)
library(dostats)
context("ID maker")
}
test_that("make_new_id", {
id_maker <- make_new_id()
expect_that(id_maker$new(), equals(1))
expect_that(id_maker$new(), equals(2))

id_maker$reset(100)
expect_that(id_maker$new(), equals(101))

id_maker$reset()
expect_that(id_maker$new(3), equals(1:3))


})
