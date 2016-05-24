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
#' S4 class Representing a Fuzzy Number with Sides Given by Power Functions
#'
#' @description
#' Bodjanova-type fuzzy numbers which sides are given by power functions
#' are defined using four coefficients
#' \code{a1} <= \code{a2} <= \code{a3} <= \code{a4},
#' and parameters \code{p.left}, \code{p.right}>0, which determine
#' exponents for the side functions.
#'
#' @details
#' We have \eqn{\mathtt{left}(x)=x^{\mathtt{p.left}}}{left(x)=x^p.left},
#' and \eqn{\mathtt{right}(x)=(1-x)^{\mathtt{p.right}}}{right(x)=(1-x)^p.right}.
#'
#' This class is a natural generalization of trapezoidal FNs.
#' For other see \linkS4class{PiecewiseLinearFuzzyNumber}.
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{a1}, \code{a2}, \code{a3}, \code{a4},
#'    \code{lower}, \code{upper}, \code{left}, \code{right}:}{
#'    Inherited from the \code{\linkS4class{FuzzyNumber}} class.}
#'    \item{\code{p.left}:}{single numeric value; 1.0 for a trapezoidal FN.}
#'    \item{\code{p.right}:}{single numeric value; 1.0 for a trapezoidal FN.}
#'  }
#'
#' @section Extends:
#' Class \code{\linkS4class{FuzzyNumber}}, directly.
#'
#' @seealso \code{\link{PowerFuzzyNumber}} for a convenient constructor,
#' \code{\link{as.PowerFuzzyNumber}} for conversion of objects to this class.
#' 
#' @references
#' Bodjanova S. (2005), Median value and median interval of a fuzzy number,
#' Information Sciences 172, pp. 73-89.
#'
#'
#' @exportClass PowerFuzzyNumber
#' @name PowerFuzzyNumber-class
#' @rdname PowerFuzzyNumber-class
#' @seealso \code{\link{PowerFuzzyNumber}} for a convenient constructor
#' @docType class
#' @family PowerFuzzyNumber-method
#' @examples
#' showClass("PowerFuzzyNumber")
#' showMethods(classes="PowerFuzzyNumber")
setClass(
   Class="PowerFuzzyNumber",
   representation(
      p.left="numeric",
      p.right="numeric"
   ),
   prototype=prototype(
      p.left=1.0,
      p.right=1.0
   ),
   validity=function(object)
   {
      if (object@p.left <= 0) return("`p.left' should be > 0")
      if (object@p.right <= 0) return("`p.right' should be > 0")

      return(TRUE);
   },
   contains="FuzzyNumber"
)



setMethod(
   f="initialize",
   signature("PowerFuzzyNumber"),
   definition=function(.Object, ...)
   {
      .Object <- callNextMethod()

      e <- new.env();

      e$p.left  <- p.left  <- .Object@p.left    # p.left <- ... to avoid
      e$p.right <- p.right <- .Object@p.right   # PKG CHECK problems

      .Object@left   <- function(x)             x^(p.left)
      .Object@right  <- function(x)         (1-x)^(p.right)
      .Object@lower  <- function(alpha)     alpha^(1.0/p.left)
      .Object@upper  <- function(alpha)   1-alpha^(1.0/p.right)

      environment(.Object@left)  <- e
      environment(.Object@right) <- e
      environment(.Object@lower) <- e
      environment(.Object@upper) <- e

      return(.Object)
   }
)


#' @title
#' Creates a Fuzzy Number with Sides Given by Power Functions
#'
#' @description
#' For convenience, objects of class \code{\linkS4class{PowerFuzzyNumber}}
#' may be created with this function.
#'
#' @param a1 a number specyfing left bound of the support
#' @param a2 a number specyfing left bound of the core
#' @param a3 a number specyfing right bound of the core
#' @param a4 a number specyfing right bound of the support
#' @param p.left a positive number specyfing the exponent for the left side
#' @param p.right a positive number specyfing the exponent for the right side
#' @return Object of class \code{\linkS4class{PowerFuzzyNumber}}
#' @export
#'
#' @family PowerFuzzyNumber-method
PowerFuzzyNumber <- function(a1, a2, a3, a4, p.left=1.0, p.right=1.0)
{
   new("PowerFuzzyNumber", a1=a1, a2=a2, a3=a3, a4=a4,
                                      p.left=p.left, p.right=p.right)
}
