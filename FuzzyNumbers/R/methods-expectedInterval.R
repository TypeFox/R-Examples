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
#' Calculate the Expected Interval of a Fuzzy Number
#'
#' @description
#' We have \eqn{EI(A) := [\int_0^1  A_L(\alpha)\,d\alpha,\int_0^1  A_U(\alpha)\,d\alpha]
#' }{EI(A) := [int_0^1  A_L(\alpha) d\alpha, int_0^1  A_U(\alpha) d\alpha]},
#' see (Duboid, Prade, 1987).
#'
#' @details
#' Note that if an instance of the \code{FuzzyNumber} or \code{DiscontinuousFuzzyNumber} class
#' is given, the calculation is performed via numerical integration.
#' Otherwise, the computation is exact.
#'
#' @param object a fuzzy number
#' @param ... for \code{FuzzyNumber} and \code{DiscontinuousFuzzyNumber} - additional arguments passed to \code{\link{integrateAlpha}}
#' 
#' @return Returns a numeric vector of length 2.
#'
#' @references
#' Dubois D., Prade H. (1987), The mean value of a fuzzy number,
#'  Fuzzy Sets and Systems 24, pp. 279-300.
#'
#' @exportMethod expectedInterval
#' @docType methods
#' @name expectedInterval
#' @family FuzzyNumber-method
#' @family TrapezoidalFuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family PowerFuzzyNumber-method
#' @rdname expectedInterval-methods
#' @aliases expectedInterval,FuzzyNumber-method
#'          expectedInterval,TrapezoidalFuzzyNumber-method
#'          expectedInterval,PiecewiseLinearFuzzyNumber-method
#'          expectedInterval,PowerFuzzyNumber-method
#' @usage
#' \S4method{expectedInterval}{FuzzyNumber}(object, ...)
#'
#' \S4method{expectedInterval}{TrapezoidalFuzzyNumber}(object)
#'
#' \S4method{expectedInterval}{PiecewiseLinearFuzzyNumber}(object)
#'
#' \S4method{expectedInterval}{PowerFuzzyNumber}(object)
setGeneric("expectedInterval",
           function(object, ...) standardGeneric("expectedInterval"))




setMethod(
   f="expectedInterval",
   signature(object="FuzzyNumber"),
   definition=function(object, ...)
   {
      if (is.na(object@lower(0))) return(c(NA_real_, NA_real_))

      return(c(
         integrateAlpha(object, "lower", 0, 1, ...),
         integrateAlpha(object, "upper", 0, 1, ...)
      ))
   }
)





setMethod(
   f="expectedInterval",
   signature(object="TrapezoidalFuzzyNumber"),
   definition=function(object)
   {
      return(0.5*c((object@a2+object@a1), (object@a4+object@a3)))
   }
)






setMethod(
   f="expectedInterval",
   signature(object="PiecewiseLinearFuzzyNumber"),
   definition=function(object)
   {
      xl <- c(object@a1, object@knot.left,  object@a2)
      xr <- c(object@a3, object@knot.right, object@a4)
      dal <- diff(c(0,     object@knot.alpha,  1))
      dar <- diff(c(1, rev(object@knot.alpha), 0))

      return(c(
         sum( (xl[-object@knot.n-2]+0.5*diff(xl))*dal ),
         sum(-(xr[-object@knot.n-2]+0.5*diff(xr))*dar )
      ))
   }
)



setMethod(
   f="expectedInterval",
   signature(object="PowerFuzzyNumber"),
   definition=function(object)
   {
      return(c(
         (object@a1+object@p.left*(object@a2-object@a1)/(object@p.left+1) ),
         (object@a3+(object@a4-object@a3)/(object@p.right+1) )
      ))
   }
)
