## This file is part of the FuzzyNumbers library.
##
## Copyright 2012-2014 Marek Gagolewski
##
##
## FuzzyNumbers is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## FuzzyNumbers is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with FuzzyNumbers. If not, see <http://www.gnu.org/licenses/>.



#' @title
#' Calculate the Width of a Fuzzy Number
#'
#' @description
#' The width (Chanas, 2001) is a measure of nonspecificity of a fuzzy number.
#'
#' @details
#' The width of \eqn{A} is defined as
#' \eqn{width(A) := EI_U(A) - EI_L(A)},
#' where \eqn{EI} is the \code{\link{expectedInterval}}.
#'
#' @param object a fuzzy number
#' @param ... additional arguments passed to \code{\link{expectedInterval}}
#' 
#' @return Returns a single numeric value.
#'
#'
#' @exportMethod width
#' @docType methods
#' @name width
#' @family FuzzyNumber-method
#' @family characteristics
#' @rdname width-methods
#' @aliases width,FuzzyNumber-method
#' @usage
#' \S4method{width}{FuzzyNumber}(object, ...)
#' @references
#' Chanas S. (2001), On the interval approximation of a fuzzy number,
#' Fuzzy Sets and Systems 122, pp. 353-356.
setGeneric("width",
           function(object, ...) standardGeneric("width"))




setMethod(
   f="width",
   signature(object="FuzzyNumber"),
   definition=function(object, ...)
   {
      return(diff(expectedInterval(object, ...)))
   }
)
