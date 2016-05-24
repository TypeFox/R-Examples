{############################################################################### 
# test-base.patterns.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 7/30/2012
# 
# DESCRIPTION
# ===========
# testing for base patterns. 
# 
# LICENSE
# ========
# lint is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# dostats is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see http://www.gnu.org/licenses/.
# 
# LOG
# ===
# 1/18/2013 Checked to work under R-3.0
# 
}###############################################################################

context("Base Patterns")

test_that("start.characters", {
    a <- c('a', 'A', '1', '.', '_')
    y <- c( T ,  T ,  F ,  T ,  F )
    z <- grepl(start.characters, a, perl=T)
    expect_equal(y, z)
})



