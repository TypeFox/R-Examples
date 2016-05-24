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
#' Calculate the Alpha-Interval of a Fuzzy Number
#'
#' @description
#' We have \eqn{\alpha-Int(A) := [\int_0^1 \alpha A_L(\alpha)\,d\alpha, \int_0^1 \alpha A_U(\alpha)\,d\alpha]
#' }{\alpha-Int(A) := [int_0^1 \alpha A_L(\alpha) d\alpha, int_0^1 \alpha A_U(\alpha) d\alpha]}.
#'
#' @details
#' Note that if an instance of the \code{FuzzyNumber} or 
#' \code{DiscontinuousFuzzyNumber} class
#' is given, the calculation is performed via numerical integration.
#' Otherwise, the computation is exact.
#'
#' @param object a fuzzy number
#' @param ... for \code{FuzzyNumber} and \code{DiscontinuousFuzzyNumber}
#'  - additional arguments passed to \code{\link{integrateAlpha}}
#' 
#' @return Returns numeric vector of length 2.
#'
#' @exportMethod alphaInterval
#' @docType methods
#' @name alphaInterval
#' @family FuzzyNumber-method
#' @family TrapezoidalFuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family PowerFuzzyNumber-method
#' @rdname alphaInterval-methods
#' @aliases alphaInterval,FuzzyNumber-method
#'          alphaInterval,TrapezoidalFuzzyNumber-method
#'          alphaInterval,PiecewiseLinearFuzzyNumber-method
#'          alphaInterval,PowerFuzzyNumber-method
#' @usage
#' \S4method{alphaInterval}{FuzzyNumber}(object, ...)
#'
#' \S4method{alphaInterval}{TrapezoidalFuzzyNumber}(object)
#'
#' \S4method{alphaInterval}{PiecewiseLinearFuzzyNumber}(object)
#'
#' \S4method{alphaInterval}{PowerFuzzyNumber}(object)
setGeneric("alphaInterval",
   function(object, ...) standardGeneric("alphaInterval"))


setMethod(
   f="alphaInterval",
   signature(object="FuzzyNumber"),
   definition=function(object, ...)
   {
      if (is.na(object@lower(0)) || is.na(object@upper(0))) return(c(NA_real_, NA_real_))

      return(c(
         integrateAlpha(object, "lower", 0, 1, weight=identity, ...),
         integrateAlpha(object, "upper", 0, 1, weight=identity, ...)
      ))
   }
)





setMethod(
   f="alphaInterval",
   signature(object="TrapezoidalFuzzyNumber"),
   definition=function(object)
   {
      return(c(
         object@a1*0.5+(object@a2-object@a1)/3,
         object@a3*0.5+(object@a4-object@a3)/6
      ))
   }
)






setMethod(
   f="alphaInterval",
   signature(object="PiecewiseLinearFuzzyNumber"),
   definition=function(object)
   {
      xl <- c(object@a1, object@knot.left,  object@a2)
      xr <- c(object@a3, object@knot.right, object@a4)
      al <- c(0,     object@knot.alpha,  1)
      ar <- c(1, rev(object@knot.alpha), 0)
      dxl <- diff(xl)
      dxr <- diff(xr)
      dal <- diff(al)
      dar <- diff(ar)

      return(c(
          sum( diff(al^2)*(xl[-object@knot.n-2]-al[-object@knot.n-2]*dxl/dal)/2+diff(al^3)*dxl/dal/3 ),
         -sum( diff(ar^2)*(xr[-object@knot.n-2]-ar[-object@knot.n-2]*dxr/dar)/2+diff(ar^3)*dxr/dar/3 )
      ))
   }
)



setMethod(
   f="alphaInterval",
   signature(object="PowerFuzzyNumber"),
   definition=function(object)
   {
      return(c(
         (2*object@a2*object@p.left+object@a1)/(4*object@p.left+2),
         (2*object@a3*object@p.right+object@a4)/(4*object@p.right+2)
      ))
   }
)
