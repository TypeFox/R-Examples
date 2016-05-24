{############################################################################### 
# styles.performance.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# pre-defined performance problems.
# 
# LICENSE
# ========
# lint is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
# 
# lint is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see http://www.gnu.org/licenses/.
# 
}###############################################################################

#' @name performance-styles
#' @rdname performance-styles
#' @docType data
#' @title Performance Enhancing Styles
#' @description
#'  This collection of styles assert checks on known performance improvements.
#' @format  Lint style checks.
#' @seealso lint
#' 
#' @include base.patterns.R
#' @exportPattern ^performance\\..*$
NULL

#' @rdname performance-styles
#' @export performance.square
performance.square <- {list(
    pattern = perl(sprintf("(%s)\\s*\\*\\s*\\1", name.pattern))
)}
.testinfo.performance.square <- {list(
    lines = c('2*2'
            , 'a*a'
            , 'a^2')
  , results = data.frame(line1=2, col1=1, line2=2, col2=3)
)}




