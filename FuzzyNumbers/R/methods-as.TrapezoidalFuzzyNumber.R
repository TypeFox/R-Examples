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
#' Converts an Object to a Trapezoidal Fuzzy Number
#'
#' @description
#' This method is only for exact conversion.
#' For other cases (e.g. general FNs), use
#' \code{\link{trapezoidalApproximation}}.
#'
#' @usage
#' \S4method{as.TrapezoidalFuzzyNumber}{numeric}(object)
#'
#' \S4method{as.TrapezoidalFuzzyNumber}{FuzzyNumber}(object)
#'
#' \S4method{as.TrapezoidalFuzzyNumber}{PowerFuzzyNumber}(object)
#'
#' \S4method{as.TrapezoidalFuzzyNumber}{PiecewiseLinearFuzzyNumber}(object)
#'
#' \S4method{as.TrapezoidalFuzzyNumber}{TrapezoidalFuzzyNumber}(object)
#'
#' @param object a fuzzy number or a single numeric value (crisp number)
#' or vector of length two (interval)
#' 
#' @return Returns an bject of class \code{\linkS4class{TrapezoidalFuzzyNumber}}.
#'
#'
#' @name as.TrapezoidalFuzzyNumber
#' @docType methods
#' @rdname as.TrapezoidalFuzzyNumber
#' @family FuzzyNumber-method
#' @family TrapezoidalFuzzyNumber-method
#' @family PowerFuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family conversion
#' @aliases as.TrapezoidalFuzzyNumber,numeric-method
#'          as.TrapezoidalFuzzyNumber,FuzzyNumber-method
#'          as.TrapezoidalFuzzyNumber,TrapezoidalFuzzyNumber-method
#'          as.TrapezoidalFuzzyNumber,PowerFuzzyNumber-method
#'          as.TrapezoidalFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' @exportMethod as.TrapezoidalFuzzyNumber
setGeneric("as.TrapezoidalFuzzyNumber",
           function(object) standardGeneric("as.TrapezoidalFuzzyNumber"))




setMethod(
   f="as.TrapezoidalFuzzyNumber",
   signature(object="TrapezoidalFuzzyNumber"),
   definition=function(object)
   {
      object # returns as-is
   })


setMethod(
   f="as.TrapezoidalFuzzyNumber",
   signature(object="FuzzyNumber"),
   definition=function(object)
   {
      stop("This method is only for exact conversion. Use trapezoidalApproximation() instead.")
   })


setMethod(
   f="as.TrapezoidalFuzzyNumber",
   signature(object="PiecewiseLinearFuzzyNumber"),
   definition=function(object)
   {
      if (object@knot.n == 0)
      {
         new("TrapezoidalFuzzyNumber",
            a1=object@a1, a2=object@a2, a3=object@a3, a4=object@a4)
      }
      else
         stop("This method is only for exact conversion. Use trapezoidalApproximation() instead.")
   })


setMethod(
   f="as.TrapezoidalFuzzyNumber",
   signature(object="PowerFuzzyNumber"),
   definition=function(object)
   {
      if ((object@p.left == 1 || object@a1 == object@a2)
          && (object@p.right == 1 || object@a3 == object@a4)) {
         new("TrapezoidalFuzzyNumber",
             a1=object@a1, a2=object@a2, a3=object@a3, a4=object@a4)
      }
      else
         stop("This method is only for exact conversion. Use trapezoidalApproximation() instead.")
   })


setMethod(
   f="as.TrapezoidalFuzzyNumber",
   signature(object="numeric"),
   definition=function(object)
   {
      stopifnot(is.finite(object))

      if (length(object) == 1)
         new("TrapezoidalFuzzyNumber",
             a1=object, a2=object, a3=object, a4=object)
      else if (length(object) == 2)
         new("TrapezoidalFuzzyNumber",
             a1=object[1], a2=object[1], a3=object[2], a4=object[2])
      else
         stop("`object` should be a numeric vector of length 1 or 2")
   })
