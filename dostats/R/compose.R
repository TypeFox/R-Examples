{###############################################################################
# compose.R
# This file is part of the R package dostats.
# 
# Copyright 2012 Andrew Redd
# Date: 6/1/2012
# 
# DESCRIPTION
# ===========
# Composition of functions.
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

#' Nest functions
#' @author Andrew Redd
#' @aliases nest compose composition
#' @param ... functions to be nested together
#' @param .list alternatively an explicit list of functions.  
#'              If specified \code{...} will be ignored.
#' @details
#' compose creates a functional composition of the listed functions.
#' Functional composition of functions f and g is defined as f(g(.)).
#' Order matters the right most function listed will be the innermost 
#' function in the composition, same with the operator version.
#' To remember the order lists will be the order read out, ie. 
#' compose(f,g) = f(g(x))
#'
#' When using the operator version it is good to remember that parentheses 
#' are recommended see the examples
#'
#' @return new function consisting of the functions nested
#' @export
#' @keywords utilities, misc
#' @examples
#' compose(any, is.na)(c(NA,1:3))
#' (sum%.%is.na)(c(1,NA))  #correct
#' \dontrun{
#' sum%.%is.an(NA)  #incorrect
#' }
compose <- function(..., .list){
  l <- if(missing(.list)) {
    list(...)
  } else {
    .list
  }
  body <- as.name('...')
  for(i in rev(l)){
    body <- as.call(list(i,body))
  }
  as.function(append(alist(...=), body))
}

#' @rdname compose
#' @usage x \%.\% y
#' @param x a function
#' @param y a function
#' @export
`%.%` <- function(x,y){
  compose(.list=list(x,y))
}


