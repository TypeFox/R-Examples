{###############################################################################
# collect.R
# This file is part of the R package dostats.
# 
# Copyright 2012 Andrew Redd
# Date: 6/1/2012
# 
# DESCRIPTION
# ===========
# collect function.  Now Deprecated, use builtin function Reduce.
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

#' collect results
#' @param v a vector, list, array, etc.
#' @param fun a function to collect on
#' @param ... passed to f
#' 
#' @details
#' Collect results by resursively calling the elements of the vector \code{v}.
#' The first two elements are called as \code{fun(v[1], v[2],...)}  The result is x.
#' Then f(x, v[3]) is called and so forth, until all elements has been exhausted.
#' 
#' as such \code{fun} must take two arguments and return a single element, although
#' there are no restrictions on what that single thing might be.
#' 
#' @export
#' @examples
#' collect(v=letters, fun=function(x,y,...)paste(y,x, ...), sep='/')
collect <- function(v, fun, ...){
  len <- length(v)
  if(len < 2)
    stop("vector v not long enough.")
  x <- fun(v[[1]], v[[2]], ...)
  if(len>2) for(i in seq(from=3, to=len, by=1))
    x <- fun(x, v[[i]], ...)
  x
}


