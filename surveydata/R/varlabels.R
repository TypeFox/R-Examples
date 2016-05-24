# Defines varlabels methods

#
#  surveydata/R/varlabels.R by Andrie de Vries  Copyright (C) 2011-2012
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


#' Returns and updates variable.labels attribute of surveydata object.
#' 
#' In a surveydata object, the \code{variable.labels} attribute store metadata about the original question text (see \code{\link[foreign]{read.spss}} for details).  The function \code{varlabels} returns the \code{variable.labels} attribute of data, and \code{varlabels(x) <- value} updates this attribute.
#' 
#' In a surveydata object, the \code{varlabels} attribute is a named character vector, where the names correspond to the names of the the columns in 
#' 
#' @param x surveydata object
#' @param value New value
#' @export  
#' @seealso \code{\link{surveydata-package}}
#' @family Attribute functions
#' @example /inst/examples/example-varlabels.R
varlabels <- function(x){
  attr(x, "variable.labels")
}

#' @rdname varlabels
#' @usage "varlabels(x) <- value"
#' @aliases varlabels<-
#' @export varlabels<-
#' @family Attribute functions
"varlabels<-" <- function(x, value){
  attr(x, "variable.labels") <- value
  x
}

