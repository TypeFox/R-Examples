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
#' Evaluate the Membership Function
#'
#' @description
#' This function returns the value(s) of the membership function
#' of a fuzzy number at given point(s).
#'
#'
#' @param object a fuzzy numbers
#' @param x numeric vector
#'
#' @return Returns a numeric vector.
#'
#' @exportMethod evaluate
#' @docType methods
#' @name evaluate
#' @family FuzzyNumber-method
#' @family membership_function
#' @rdname evaluate-methods
#' @aliases evaluate,FuzzyNumber,numeric-method
#' @usage
#' \S4method{evaluate}{FuzzyNumber,numeric}(object, x)
#' @examples
#' T <- TrapezoidalFuzzyNumber(1,2,3,4)
#' evaluate(T, seq(0, 5, by=0.5))
setGeneric("evaluate",
           function(object, x) standardGeneric("evaluate"))



setMethod(
   f="evaluate",
   signature(object="FuzzyNumber", x="numeric"),
   definition=function(object, x)
   {
      y <- numeric(length(x)) # prealloc
      wh1 <- which(x >= object@a1 & x <  object@a2)
      y[wh1] <- object@left ((x[wh1]-object@a1)/(object@a2-object@a1))
      wh2 <- which(x >  object@a3 & x <= object@a4)
      y[wh2] <- object@right((x[wh2]-object@a3)/(object@a4-object@a3))
      y[x >= object@a2 & x <= object@a3] <- 1.0
      names(y) <- format(x)
      y
   }
)
