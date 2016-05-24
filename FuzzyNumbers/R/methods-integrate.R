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


# internal - test
integrate_identity <- function(object, which=c("lower", "upper"), from, to, n=100L)
{
   which <- match.arg(which)
   stopifnot(is.numeric(from), length(from) == 1, is.finite(from))
   stopifnot(is.numeric(to), length(to) == 1, is.finite(to))
   stopifnot(is.numeric(n), length(n) == 1, n > 1)

   # Boole's rule (Newton-Cotes of degree 4)
   br <- seq(from, to, length.out=4*n-3)
   if (which == 'lower')
      fval <- (object@a1+(object@a2-object@a1)*object@lower(br))*br
   else
      fval <- (object@a3+(object@a4-object@a3)*object@upper(br))*br

   sum(fval * c(7, rep(c(32, 12, 32, 14), times=n-2), 32, 12, 32, 7)/90)*(to-from)/(n-1)
}


#' @title
#' Integrate a Function with at Most Finite Number of Discontinuities *EXPERIMENTAL*
#'
#' @description
#' The function uses multiple calls to \code{\link{integrate}}.
#'
#' @param f an R function taking a numeric vector of length 1 as its first
#'         argument and returning a numeric vector of length 1
#' @param from the lower limit of integration
#' @param to the upper limit of integration
#' @param discontinuities nondecreasingly sorted numeric vector which indicates
#'          the points at which \code{f} is discontinuous
#' @param ... further arguments to be passed to the \code{\link{integrate}} function.
#' 
#' @return Returns the estimate of the integral.
#' @export
integrate_discont_val <- function(f, from, to, discontinuities=numeric(0), ...)
{
   if (!is.numeric(discontinuities))
      stop("`discontinuities' should be numeric")

   stopifnot(from <= to)

   discontinuities <- discontinuities[discontinuities > from & discontinuities < to]
   m <- length(discontinuities)

   if (m == 0)
      return(integrate(f=f, lower=from, upper=to, ...)$value)

   if (is.unsorted(discontinuities))
      stop("`discontinuities' should be ordered nondecreasingly")

   x <- c(from, discontinuities, to)
   v <- numeric(m+1)

   for (i in 1:(m+1))
   {
      v[i] <- integrate(f=f, lower=x[i]+.Machine$double.eps^0.5, upper=x[i+1]-.Machine$double.eps^0.5, ...)$value
   }

   return(sum(v))
}







#' @title
#' Numerically Integrate Alpha-Cut Bounds
#'
#' @description
#' Integrates numerically a transformed or weighted lower or upper alpha-cut bound of a fuzzy number.
#'
#' @param object a fuzzy number
#' @param which one of \code{"lower"}, \code{"upper"}
#' @param from numeric
#' @param to numeric
#' @param weight a function or NULL
#' @param transform a function or NULL
#' @param ... additional arguments passed to \code{\link{integrate}} or \code{\link{integrate_discont_val}}
#' 
#' @return Returns a single numeric value.
#'
#' @exportMethod integrateAlpha
#' @docType methods
#' @name integrateAlpha
#' @family FuzzyNumber-method
#' @family DiscontinuousFuzzyNumber-method
#' @rdname integrateAlpha-methods
#' @aliases integrateAlpha,FuzzyNumber,character,numeric,numeric-method
#'          integrateAlpha,DiscontinuousFuzzyNumber,character,numeric,numeric-method
#' @usage
#' \S4method{integrateAlpha}{FuzzyNumber,character,numeric,numeric}(object, which=c("lower", "upper"),
#'    from=0, to=1, weight=NULL, transform=NULL, ...)
#'
#' \S4method{integrateAlpha}{DiscontinuousFuzzyNumber,character,numeric,numeric}(object, which=c("lower", "upper"),
#'    from=0, to=1, weight=NULL, transform=NULL, ...)
setGeneric("integrateAlpha",
     function(object, which, from, to, ...) standardGeneric("integrateAlpha"))




setMethod(
   f="integrateAlpha",
   signature(object="FuzzyNumber", which="character",
             from="numeric",       to="numeric"),
   definition=function(object, which=c("lower","upper"),
                       from=0, to=1, weight=NULL, transform=NULL, ...)
   {
      which <- match.arg(which)

      if (length(from) != 1 || length(to) != 1 || from < 0 || to > 1)
         stop("invalid `from' or `to'")

      if (!is.null(weight) && (class(weight) != "function" || length(formals(weight)) != 1))
         stop("`weight' should be a function with 1 parameter")

      if (!is.null(transform) && (class(transform) != "function" || length(formals(transform)) != 2))
         stop("`transform' should be a function with 2 parameter")

      if (!is.null(weight) && !is.null(transform))
         stop("specify either `weight', `transform' or none")

      if (which == "lower")
      {
         if (!is.null(weight))
         {
            fun <- function(alpha)
               (object@a1+(object@a2-object@a1)*object@lower(alpha))*weight(alpha)
         } else if (!is.null(transform))
         {
            fun <- function(alpha)
               transform(alpha, object@a1+(object@a2-object@a1)*object@lower(alpha))
         } else
         {
            fun <- function(alpha)
               object@a1+(object@a2-object@a1)*object@lower(alpha)
         }
      } else
      {
         if (!is.null(weight))
         {
            fun <- function(alpha)
               (object@a3+(object@a4-object@a3)*object@upper(alpha))*weight(alpha)
         } else if (!is.null(transform))
         {
            fun <- function(alpha)
               transform(alpha, object@a3+(object@a4-object@a3)*object@upper(alpha))
         } else
         {
            fun <- function(alpha)
               object@a3+(object@a4-object@a3)*object@upper(alpha)
         }
      }

      integrate(f=fun, from, to, ...)$value
   }
)


setMethod(
   f="integrateAlpha",
   signature(object="DiscontinuousFuzzyNumber", which="character",
             from="numeric",       to="numeric"),
   definition=function(object, which=c("lower", "upper"),
      from=0, to=1, weight=NULL, transform=NULL, ...)
   {
      which <- match.arg(which)

      if (length(from) != 1 || length(to) != 1 || from < 0 || to > 1)
         stop("invalid `from' or `to'")

      if (!is.null(weight) && (class(weight) != "function" || length(formals(weight)) != 1))
         stop("`weight' should be a function with 1 parameter")

      if (!is.null(transform) && (class(transform) != "function" || length(formals(transform)) != 2))
         stop("`transform' should be a function with 2 parameter")

      if (!is.null(weight) && !is.null(transform))
         stop("specify either `weight', `transform' or none")

      if (which == "lower")
      {
         if (!is.null(weight))
         {
            fun <- function(alpha)
               (object@a1+(object@a2-object@a1)*object@lower(alpha))*weight(alpha)
         } else if (!is.null(transform))
         {
            fun <- function(alpha)
               transform(alpha, object@a1+(object@a2-object@a1)*object@lower(alpha))
         } else
         {
            fun <- function(alpha)
               object@a1+(object@a2-object@a1)*object@lower(alpha)
         }

         disconts <- object@discontinuities.lower

      } else
      {
         if (!is.null(weight))
         {
            fun <- function(alpha)
               (object@a3+(object@a4-object@a3)*object@upper(alpha))*weight(alpha)
         } else if (!is.null(transform))
         {
            fun <- function(alpha)
               transform(alpha, object@a3+(object@a4-object@a3)*object@upper(alpha))
         } else
         {
            fun <- function(alpha)
               object@a3+(object@a4-object@a3)*object@upper(alpha)
         }

         disconts <- object@discontinuities.upper
      }

      integrate_discont_val(fun, from, to,
         discontinuities=disconts, ...)
   }
)
