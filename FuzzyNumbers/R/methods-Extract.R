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
#' FuzzyNumber Slot Accessors
#'
#' @description
#' For possible slot names see man pages for the \linkS4class{FuzzyNumber} class and its derivatives
#'
#' @details
#' All slot accessors are read-only.
#'
#' @param x a fuzzy number
#' @param i character; slot name
#'
#' @return Returns the slot value.
#'
#' @usage
#' \S4method{[}{FuzzyNumber,character}(x, i)
#'
#' \S4method{[}{PiecewiseLinearFuzzyNumber,character}(x, i)
#'
#' \S4method{[}{PowerFuzzyNumber,character}(x, i)
#'
#' \S4method{[}{DiscontinuousFuzzyNumber,character}(x, i)
#'
#' @name Extract
#' @rdname Extract-methods
#' @docType methods
#' @family FuzzyNumber,character-method
#' @family PiecewiseLinearFuzzyNumber,character-method
#' @family PowerFuzzyNumber,character-method
#' @family DiscontinuousFuzzyNumber,character-method
#' @exportMethod [
#' @aliases [,FuzzyNumber,character-method
#'          [,PiecewiseLinearFuzzyNumber,character-method
#'          [,PowerFuzzyNumber,character-method
#'          [,DiscontinuousFuzzyNumber,character-method
#' @examples
#' A <- FuzzyNumber(1,2,3,4)
#' A["a1"]
#' A["right"]
invisible(NULL)



setMethod(
   f="[",
   signature(x="FuzzyNumber", i="character"),
   definition=function(x, i)
   {
      switch(i,
             "a1"    = x@a1,
             "a2"    = x@a2,
             "a3"    = x@a3,
             "a4"    = x@a4,
             "left"  = x@left,
             "right" = x@right,
             "lower" = x@lower,
             "upper" = x@upper,
             NULL
      )
   }
)


setMethod(
   f="[",
   signature(x="PiecewiseLinearFuzzyNumber", i="character"),
   definition=function(x, i)
   {
      switch(i,
             "a1"    = x@a1,
             "a2"    = x@a2,
             "a3"    = x@a3,
             "a4"    = x@a4,
             "left"  = x@left,
             "right" = x@right,
             "lower" = x@lower,
             "upper" = x@upper,
             "knot.n"=x@knot.n,
             "knot.alpha"=x@knot.alpha,
             "knot.left"=x@knot.left,
             "knot.right"=x@knot.right,
             "knots"=matrix(c(x@knot.alpha, x@knot.left, rev(x@knot.right)),
                            ncol=3,
                            dimnames=list(
                               paste("knot_", 1:x@knot.n, sep=""),
                               c("alpha", "L", "U")
                            )
             ),
             "allknots"=matrix(c(0,x@knot.alpha,1,  x@a1, x@knot.left, x@a2, x@a4, rev(x@knot.right), x@a3),
                               ncol=3,
                               dimnames=list(
                                  c("supp", paste("knot_", 1:x@knot.n, sep=""), "core"),
                                  c("alpha", "L", "U")
                               )
             ),
             NULL
      )
      #     return(callNextMethod()) # does not work...
   }
)


setMethod(
   f="[",
   signature(x="PowerFuzzyNumber", i="character"),
   definition=function(x, i)
   {
      switch(i,
             "a1"    = x@a1,
             "a2"    = x@a2,
             "a3"    = x@a3,
             "a4"    = x@a4,
             "left"  = x@left,
             "right" = x@right,
             "lower" = x@lower,
             "upper" = x@upper,
             "p.left"  = x@p.left,
             "p.right" = x@p.right,
             NULL
      )
   }
)




setMethod(
   f="[",
   signature(x="DiscontinuousFuzzyNumber", i="character"),
   definition=function(x, i)
   {
      switch(i,
             "a1"    = x@a1,
             "a2"    = x@a2,
             "a3"    = x@a3,
             "a4"    = x@a4,
             "left"  = x@left,
             "right" = x@right,
             "lower" = x@lower,
             "upper" = x@upper,
             "discontinuities.left"  = x@discontinuities.left,
             "discontinuities.right" = x@discontinuities.right,
             "discontinuities.lower" = x@discontinuities.lower,
             "discontinuities.upper" = x@discontinuities.upper,
             NULL
      )
   }
)
