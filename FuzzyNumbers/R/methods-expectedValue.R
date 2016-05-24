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
#' Calculate the Expected Value of a Fuzzy Number
#'
#' @description
#' The calculation of the so-called expected value is one of possible methods to
#' deffuzify a fuzzy number.
#'
#' @details
#' The expected value of \eqn{A} is defined as
#' \eqn{EV(A) := (EI_U(A) + EI_L(A))/2}{EV(A) := (EI_U(A) + EI_L(A))/2},
#' where \eqn{EI} is the \code{\link{expectedInterval}.}
#'
#' @param object a fuzzy number
#' @param ... additional arguments passed to \code{\link{expectedInterval}}
#' 
#' @return Returns a single numeric value.
#'
#' @exportMethod expectedValue
#' @docType methods
#' @name expectedValue
#' @family FuzzyNumber-method
#' @family deffuzification
#' @family characteristics
#' @rdname expectedValue-methods
#' @aliases expectedValue,FuzzyNumber-method
#' @usage
#' \S4method{expectedValue}{FuzzyNumber}(object, ...)
setGeneric("expectedValue",
           function(object, ...) standardGeneric("expectedValue"))



setMethod(
   f="expectedValue",
   signature(object="FuzzyNumber"),
   definition=function(object, ...)
   {
      return(mean(expectedInterval(object, ...)))
   }
)
