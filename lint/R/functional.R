{############################################################################### 
# functional.R
# This file is part of the R package lint.
# 
# Copyright 2012 Andrew Redd
# Date: 6/16/2012
# 
# DESCRIPTION
# ===========
# functional formated tests
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


#' @name functional.Rd
#' @title function based tests
#' @aliases check_functional function_style
#'
#' @description
#'  Function based tests are much more flexible than pattern based tests and as 
#'  such more responsibility is placed on the authors of such tests.
#'  Theoretically, function based tests are capable of any type of style tests.
#'  
#'  Inclusion/Exclusion is handled by the calling function(s )and does not need
#'  to be included in the function but may be
#'  unless more complicated than include.regions
#'  and exclude.regions can handle.
#'  
#' @param f function that performs the test, the remaining parameters must be 
#'            accepted by f.
#' @param file file name of the file being tested.  File is not guaranteed to 
#'             be an on disk file but could be a \link{connection}.
#' @param lines character vector of the lines of the file.
#' @param parse.data data frame with parse data from \code{\link{getParseData}}.
#' @param ... discarded.
#' 
#' @return The return value of f shouls be a find formated data frame.
#'         in the case that the files pass the test the function should
#'         return an \link{empty.find}.
#'
#' @export check_functional
check_functional <- function(f, ..., file, lines, parse.data){
    result <- f(file=file, lines=lines, parse.data=parse.data)
    if(isTRUE(result))
        result <- empty.find
    stopifnot(valid_find(result))
    return(result)
}

