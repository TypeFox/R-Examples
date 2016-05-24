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
#' Converts an Object to a Fuzzy Number
#'
#' @description
#' Please note that applying this function on a \linkS4class{FuzzyNumber}
#' child class causes information loss, as it drops all additional slots
#' defined in the child classes.
#' \code{\linkS4class{FuzzyNumber}} is the base class for all FNs.
#' Note that some functions for TFNs or PLFNs
#' work much faster and are more precise. This function shouldn't be
#' used in normal computations.
#'
#' @usage
#' \S4method{as.FuzzyNumber}{numeric}(object)
#'
#' \S4method{as.FuzzyNumber}{FuzzyNumber}(object)
#'
#' @param object a fuzzy number or a single numeric value (crisp number)
#' or vector of length two (interval)
#' 
#' @return Returns an bject of class \code{\linkS4class{FuzzyNumber}}.
#'
#'
#' @name as.FuzzyNumber
#' @docType methods
#' @rdname as.FuzzyNumber
#' @family FuzzyNumber-method
#' @family conversion
#' @aliases as.FuzzyNumber,numeric-method
#'          as.FuzzyNumber,FuzzyNumber-method
#' @exportMethod as.FuzzyNumber
setGeneric("as.FuzzyNumber",
           function(object) standardGeneric("as.FuzzyNumber"))




setMethod(
   f="as.FuzzyNumber",
   signature(object="FuzzyNumber"),
   definition=function(object)
   {
      new("FuzzyNumber",
          a1=object@a1, a2=object@a2, a3=object@a3, a4=object@a4,
          left=object@left, right=object@right,
          lower=object@lower, upper=object@upper)
   })



setMethod(
   f="as.FuzzyNumber",
   signature(object="numeric"),
   definition=function(object)
   {
      stopifnot(is.finite(object))

      if (length(object) == 1)
         new("FuzzyNumber",
             a1=object, a2=object, a3=object, a4=object,
             left=function(x) x, right=function(x) 1-x,
             lower=function(x) x, upper=function(x) 1-x)
      else if (length(object) == 2)
         new("FuzzyNumber",
             a1=object[1], a2=object[1], a3=object[2], a4=object[2],
             left=function(x) x, right=function(x) 1-x,
             lower=function(x) x, upper=function(x) 1-x)
      else
         stop("`object` should be a numeric vector of length 1 or 2")
   })
