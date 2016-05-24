{###############################################################################
# dostats.R
# This file is part of the R package dostats
# 
# Copyright 2012 Andrew Redd
# Date: 5/30/2012
# 
# DESCRIPTION
# ===========
# dostats is a helper function for computing descriptive tables.
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

#' Convenient interface for computing statistics on a vector
#'  @author Andrew Redd
#'
#'  @param x the vector
#'  @param ... statistics to compute, must take a vector and return a vector
#'  @param .na.action the action to take on NA values, for all statistics
#'
#'  @return A one row \code{data.frame} with columns named as in \code{...}
#'  @export
#'  @seealso \code{\link[plyr]{ldply}}
#'  @keywords utilities, misc
#'  @example inst/ex_dostats.R
dostats <- function(x, ..., .na.action=na.fail){
  if(any(is.na(x)))
    x <- .na.action(x)
  funs   <- list(...)
  fnames <- names(funs)
  inames <- as.character(substitute(c(...)))[-1]
  fnames <- if(is.null(fnames)) inames else ifelse(fnames != "", fnames, inames)
  l <- structure(lapply(funs, do.call, list(x)), names=fnames)
  l <- lapply(l, function(y)if(length(y)==1) y else t(y))
  do.call(data.frame, l)
}

#' Filter by class
#' @param x vector of any class
#' @param .class string for class to filter by
#' @param ... passed to \code{\link{dostats}}
#' @return data frame of computed statistics if x is of class \code{.class}
#'         otherwise returns \code{NULL}.
#'  @export
#' @seealso \code{\link{dostats}}
class.stats <- function(.class){
  if(class(.class)!="character")
    .class=as.character(substitute(.class))
  function(x, ...){if(inherits(x, .class))
    dostats(x, ...)
  else NULL
  }
}
#' @rdname class.stats
#'  @export
numeric.stats <- class.stats(numeric)
#' @rdname class.stats
#'  @export
factor.stats  <- class.stats(factor)
#' @rdname class.stats
#'  @export
integer.stats <- class.stats(integer)