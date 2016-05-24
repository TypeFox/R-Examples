
#
#  surveydata/R/pattern.R by Andrie de Vries  Copyright (C) 2011-2012
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



#' Returns and updates pattern attribute.
#' 
#' The pattern attribute contains information about the separator character used to name subquestions in the data.  Survey software typically makes use of underscores to distinguish subquestions in a grid of questions, e.g. Q4_1, Q4_2, Q4_3, Q4_other. The function \code{pattern} returns the pattern attribute, and \code{pattern<-} updates the attribute.
#' 
#' 
#' @aliases pattern pattern<-
#' @param x surveydata object
#' @export pattern 
#' @family Attribute functions
#' @seealso \code{\link{as.surveydata}}, \code{\link{which.q}}
#' @example inst/examples/example-pattern.R
pattern <- function(x){
  attr(x, "pattern")
}

#' @rdname pattern
#' @usage pattern(x) <- value
#' @param value New value
#' @export pattern<-
"pattern<-" <- function(x, value){
  attr(x, "pattern") <- value
  x
}


#' Removes pattern from attributes list.
#'
#' @param x Surveydata object
#' @keywords Internal
rm.pattern <- function(x){
  pattern(x) <- NULL
  x
}

#' Removes pattern and variable.labels from attributes list.
#'
#' @param x Surveydata object
#' @keywords Internal
rm.attrs <- function(x){
  attr(x, "pattern") <- NULL
  attr(x, "variable.labels") <- NULL
  x
}

