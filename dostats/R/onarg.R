{###############################################################################
# onarg.R
# This file is part of the R package dostats
# 
# Copyright 2012 Andrew Redd
# Date: 6/1/2012
# 
# DESCRIPTION
# ===========
# Convenience functions for manipulating default arguments of functions.
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

#' change first argument of a function
#' @param f the function
#' @param arg the arg to be called as the first argument
#' 
#' @return a function that calls \code{f} with \code{arg} as the first argument.
#' @seealso \code{\link{wargs}},  \code{\link{dostats}}, and \code{\link{apply}}
#' @export
#' @examples
#' formals(runif)
#' onarg(runif, 'max')(1:10, 1)
#' onarg(runif, 'max')(1:10, 10)
#' #another version of contains
#' onarg(`%in%`, 'table')(letters, 'y')
onarg <- function(f, arg){
  carg <- as.character(substitute(arg))
  function(x, ...){
    do.call(f,args=append(structure(list(x), names=arg), list(...)))
  }
}

#' Does a table contain a value
#' @rdname contains
#' @aliases contains
#' @param table a table of values
#' @param y a value
#' @usage table \%contains\% y
#' @usage contains(table,y)
#' 
#' @details
#' Literally %in% in reverse order, just for convenience.
#' 
#' @return a logical vector of the same length as \code{y} indicating if 
#'   \code{y} is in \code{table}, i.e. the \code{table} contains \code{y}.
#' @seealso \code{\link{match}}
#' @export contains
#' @export
`%contains%` <- contains <- function(table,y){y %in% table}
