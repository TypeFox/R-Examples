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
#' Converts an Object to a Power Fuzzy Number
#'
#' @description
#' This method is only for exact conversion.
#'
#' @usage
#' \S4method{as.PowerFuzzyNumber}{numeric}(object)
#'
#' \S4method{as.PowerFuzzyNumber}{FuzzyNumber}(object)
#'
#' \S4method{as.PowerFuzzyNumber}{PowerFuzzyNumber}(object)
#'
#' \S4method{as.PowerFuzzyNumber}{PiecewiseLinearFuzzyNumber}(object)
#'
#' \S4method{as.PowerFuzzyNumber}{TrapezoidalFuzzyNumber}(object)
#'
#' @param object a fuzzy number or a single numeric value (crisp number)
#' or vector of length two (interval)
#' 
#' @return Returns an object of class \code{\linkS4class{PowerFuzzyNumber}}.
#'
#'
#' @name as.PowerFuzzyNumber
#' @docType methods
#' @rdname as.PowerFuzzyNumber
#' @family FuzzyNumber-method
#' @family TrapezoidalFuzzyNumber-method
#' @family PowerFuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family conversion
#' @aliases as.PowerFuzzyNumber,numeric-method
#'          as.PowerFuzzyNumber,FuzzyNumber-method
#'          as.PowerFuzzyNumber,TrapezoidalFuzzyNumber-method
#'          as.PowerFuzzyNumber,PowerFuzzyNumber-method
#'          as.PowerFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' @exportMethod as.PowerFuzzyNumber
setGeneric("as.PowerFuzzyNumber",
           function(object) standardGeneric("as.PowerFuzzyNumber"))




setMethod(
   f="as.PowerFuzzyNumber",
   signature(object="PowerFuzzyNumber"),
   definition=function(object)
   {
      object # returns as-is
   })


setMethod(
   f="as.PowerFuzzyNumber",
   signature(object="FuzzyNumber"),
   definition=function(object)
   {
      stop("This method is only for exact conversion")
   })


setMethod(
   f="as.PowerFuzzyNumber",
   signature(object="PiecewiseLinearFuzzyNumber"),
   definition=function(object)
   {
      if (object@knot.n == 0)
      {
         new("PowerFuzzyNumber",
            a1=object@a1, a2=object@a2, a3=object@a3, a4=object@a4, p.left=1, p.right=1)
      }
      else
         stop("This method is only for exact conversion")
   })


setMethod(
   f="as.PowerFuzzyNumber",
   signature(object="TrapezoidalFuzzyNumber"),
   definition=function(object)
   {
      new("PowerFuzzyNumber",
          a1=object@a1, a2=object@a2, a3=object@a3, a4=object@a4, p.left=1, p.right=1)
   })


setMethod(
   f="as.PowerFuzzyNumber",
   signature(object="numeric"),
   definition=function(object)
   {
      stopifnot(is.finite(object))

      if (length(object) == 1)
         new("PowerFuzzyNumber",
             a1=object, a2=object, a3=object, a4=object, p.left=1, p.right=1)
      else if (length(object) == 2)
         new("PowerFuzzyNumber",
             a1=object[1], a2=object[1], a3=object[2], a4=object[2], p.left=1, p.right=1)
      else
         stop("`object` should be a numeric vector of length 1 or 2")
   })
