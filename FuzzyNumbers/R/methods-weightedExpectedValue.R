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
#' Calculate the Weighted Expected Value of a Fuzzy Number
#'
#' @description
#' The calculation of the so-called weighted expected value is one of possible methods to
#' deffuzify a fuzzy number.
#'
#' For \eqn{w=0.5} we get the ordinary \code{\link{expectedValue}}.
#'
#' @details
#' The weighted expected value of \eqn{A} is defined as
#' \eqn{EV_w(A) := (1-w) EI_L(A) + w EI_U(A)}{EV_w(A) := (1-w) EI_L(A) + w EI_U(A)},
#' where \eqn{EI} is the \code{\link{expectedInterval}.}
#'
#' @param object a fuzzy number
#' @param w a single numeric value in [0,1]
#' @param ... additional arguments passed to \code{\link{expectedInterval}}
#' 
#' @return Returns a single numeric value.
#'
#' @exportMethod weightedExpectedValue
#' @docType methods
#' @name weightedExpectedValue
#' @family FuzzyNumber-method
#' @family deffuzification
#' @family characteristics
#' @rdname weightedExpectedValue-methods
#' @aliases weightedExpectedValue,FuzzyNumber-method
#' @usage
#' \S4method{weightedExpectedValue}{FuzzyNumber}(object, w=0.5, ...)
setGeneric("weightedExpectedValue",
           function(object, ...) standardGeneric("weightedExpectedValue"))




setMethod(
   f="weightedExpectedValue",
   signature(object="FuzzyNumber"),
   definition=function(object, w=0.5, ...)
   {
      stopifnot(is.numeric(w), length(w) == 1, w >= 0, w <= 1)
      EI <- expectedInterval(object, ...)
      return((1-w)*EI[1] + w*EI[2])
   }
)
