### This file is part of 'PGRdup' package for R.

### Copyright (C) 2014, ICAR-NBPGR.
#
# PGRdup is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# PGRdup is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Prints summary of \code{KWIC} object.
#' 
#' \code{print.KWIC} prints to console the summary of an object of class 
#' \code{KWIC} including the database fields(columns) used, the total number of
#' keywords and the number of distinct keywords in the index.
#' 
#' @param x An object of class \code{KWIC}.
#' @param ... Unused
#' @seealso \code{\link[PGRdup]{KWIC}}
#'   
#' @export
print.KWIC <- function(x,...) {
  cat(paste("KWIC fields : ", sep = ""))
  cat(paste(x$Fields, sep = ""))
  cat(paste("\n", "Number of keywords : ", nrow(x[[1]]), sep = ""))
  cat(paste("\n", "Number of distinct keywords : ",
            length(unique(x[[1]]$KEYWORD)), sep = ""))
}
