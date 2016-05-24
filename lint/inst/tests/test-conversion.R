{############################################################################### 
# test-base.patterns.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 7/30/2012
# 
# DESCRIPTION
# ===========
# testing for conversion between formats used in lint. 
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

context("Conversion")
source(find_example("ex-conversion.R", package='lint'))

test_that('locate2find converts to find',{
    expected <- data.frame(line1=4, col1=2, line2=4, col2=6, row.names=1)
    expect_that(
        locate2find(l)
        , is_equivalent_to(expected)
    )
})
test_that('locate2find produces valid find data.',{
    expect_that(valid_find(locate2find(l)), is_true())
})
test_that('merge_find with empty arguments returns empty find.', {
    expect_identical(merge_find(), empty.find)
})
test_that('merge_find with single argument return single argument',{
    expect_that(merge_find(f), is_identical_to(f))
})
test_that('merge_find handles many empty finds.',{
    expect_that(merge_find(empty.find, empty.find), is_identical_to(empty.find))
})
test_that('parse2find collapses data frames.', {
    in.pd <- data.frame( line1 = c(1L, 2L), col1 = c(1L, 2L)
                       , line2 = c(1L, 3L), col2 = c(1L, 3L))
    out.f <- data.frame( line1 =   1L     , col1 =   1L
                       , line2 =       3L , col2 =       3L)
    expect_that(parse2find(in.pd), equals(out.f))    
})

