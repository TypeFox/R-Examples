#  Copyright 2013, 2014 Christian Sigg
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' Cardinality of Column Vectors
#' 
#' Computes the cardinality (the sum of non-zero elements) of each column of 
#' the matrix \eqn{\mathbf{w}}{w}.
#' 
#' @export
#' @param w a numeric matrix, e.g. \code{xcoef} as returned by \code{\link{nscancor}}
cardinality <- function(w) {
    return(colSums(abs(as.matrix(w)) > 0))
}