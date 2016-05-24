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
#' Converts an Object to a Piecewise Linear Fuzzy Number
#'
#' @description
#' This method is only for exact conversion.
#' For other cases (e.g. general FNs), use
#' \code{\link{piecewiseLinearApproximation}}.
#'
#' @usage
#' \S4method{as.PiecewiseLinearFuzzyNumber}{TrapezoidalFuzzyNumber}(object, knot.n=0,
#'    knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
#'
#' \S4method{as.PiecewiseLinearFuzzyNumber}{numeric}(object, knot.n=0,
#'    knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
#'
#' \S4method{as.PiecewiseLinearFuzzyNumber}{FuzzyNumber}(object, knot.n=0,
#'    knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
#'
#' \S4method{as.PiecewiseLinearFuzzyNumber}{PiecewiseLinearFuzzyNumber}(object, knot.n=0,
#'    knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
#'
#' @param object a fuzzy number or a single numeric value (crisp number)
#' or vector of length two (interval)
#' @param knot.n the number of knots
#' @param knot.alpha \code{knot.n} alpha-cut values at knots,
#' defaults to uniformly distributed knots
#' 
#' @return Returns an bject of class \code{\linkS4class{PiecewiseLinearFuzzyNumber}}.
#'
#'
#' @name as.PiecewiseLinearFuzzyNumber
#' @docType methods
#' @rdname as.PiecewiseLinearFuzzyNumber
#' @family TrapezoidalFuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family FuzzyNumber-method
#' @family conversion
#' @aliases as.PiecewiseLinearFuzzyNumber,TrapezoidalFuzzyNumber-method
#'          as.PiecewiseLinearFuzzyNumber,numeric-method
#'          as.PiecewiseLinearFuzzyNumber,FuzzyNumber-method
#'          as.PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' @exportMethod as.PiecewiseLinearFuzzyNumber
setGeneric("as.PiecewiseLinearFuzzyNumber",
           function(object, ...) standardGeneric("as.PiecewiseLinearFuzzyNumber"))




setMethod(
   f="as.PiecewiseLinearFuzzyNumber",
   signature(object="TrapezoidalFuzzyNumber"),
   definition=function(object, knot.n=0, knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
   {
      stopifnot(length(knot.n) == 1, knot.n >= 0)
      stopifnot(knot.n == length(knot.alpha))

      a <- alphacut(object, c(0, knot.alpha, 1))
      PiecewiseLinearFuzzyNumber(knot.left=a[,1], knot.right=rev(a[,2]), knot.alpha=knot.alpha)
   })



setMethod(
   f="as.PiecewiseLinearFuzzyNumber",
   signature(object="numeric"),
   definition=function(object, knot.n=0, knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
   {
      stopifnot(is.finite(object))
      stopifnot(length(knot.n) == 1, knot.n >= 0)
      stopifnot(knot.n == length(knot.alpha))

      if (length(object) == 1)
         PiecewiseLinearFuzzyNumber(a1=object, a2=object, a3=object, a4=object,
                                    knot.n=knot.n, knot.alpha=knot.alpha,
                                    knot.left=rep(object, knot.n),
                                    knot.right=rep(object, knot.n))
      else if (length(object) == 2)
         PiecewiseLinearFuzzyNumber(a1=object[1], a2=object[1], a3=object[2], a4=object[2],
                                    knot.n=knot.n, knot.alpha=knot.alpha,
                                    knot.left=rep(object, knot.n),
                                    knot.right=rep(object, knot.n))
      else
         stop("`object` should be a numeric vector of length 1 or 2")
   })



setMethod(
   f="as.PiecewiseLinearFuzzyNumber",
   signature(object="FuzzyNumber"),
   definition=function(object, knot.n=0, knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
   {
      stop("This method is only for exact conversion. Use piecewiseLinearApproximation() instead.")
   })



setMethod(
   f="as.PiecewiseLinearFuzzyNumber",
   signature(object="PiecewiseLinearFuzzyNumber"),
   definition=function(object, knot.n=0, knot.alpha=seq(0, 1, length.out=knot.n+2)[-c(1,knot.n+2)])
   {
      stopifnot(length(knot.n) == 1, knot.n >= 0)
      stopifnot(knot.n == length(knot.alpha))

      if (length(setdiff(object@knot.alpha, knot.alpha)) != 0)
         stop("This method is only for exact conversion. Use piecewiseLinearApproximation() instead.")

      a <- alphacut(object, c(0, knot.alpha, 1))
      PiecewiseLinearFuzzyNumber(knot.left=a[,1], knot.right=rev(a[,2]))
   })
