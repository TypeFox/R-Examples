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
#' S4 class Representing a Trapezoidal Fuzzy Number
#'
#' @description
#' Trapezoidal Fuzzy Numbers have linear side functions and alpha-cut bounds.
#'
#' @details
#' Trapezoidal fuzzy numbers are among the simplest FNs.
#' Despite their simplicity, however, they include triangular FNs,
#' ``crisp'' real intervals, and ``crisp'' reals.
#' Please note that currently no separate classes for these particular TFNs types
#' are implemented in the package.
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{a1}, \code{a2}, \code{a3}, \code{a4},
#'    \code{lower}, \code{upper}, \code{left}, \code{right}:}{
#'    Inherited from the \code{\linkS4class{FuzzyNumber}} class.}
#'  }
#'
#' @section Extends:
#' Class \code{\linkS4class{FuzzyNumber}}, directly.
#'
#' @seealso \code{\link{TrapezoidalFuzzyNumber}} for a convenient constructor,
#' \code{\link{as.TrapezoidalFuzzyNumber}} for conversion of objects to this class,
#' and \code{\link{trapezoidalApproximation}} for approximation routines.
#'
#' @exportClass TrapezoidalFuzzyNumber
#' @name TrapezoidalFuzzyNumber-class
#' @rdname TrapezoidalFuzzyNumber-class
#' @docType class
#' @family TrapezoidalFuzzyNumber-method
#' @examples
#' showClass("TrapezoidalFuzzyNumber")
#' showMethods(classes="TrapezoidalFuzzyNumber")
setClass(
   Class="TrapezoidalFuzzyNumber",
   prototype=prototype(
      left=function(x) x,
      right=function(x) 1-x,
      lower=function(alpha) alpha,
      upper=function(alpha) 1-alpha
   ),
   contains="FuzzyNumber"
)


#' @title
#' Creates a Trapezoidal Fuzzy Number
#'
#' @description
#' For convenience, objects of class \code{\linkS4class{TrapezoidalFuzzyNumber}}
#' may be created with this function.
#'
#' @param a1 a number specyfing left bound of the support
#' @param a2 a number specyfing left bound of the core
#' @param a3 a number specyfing right bound of the core
#' @param a4 a number specyfing right bound of the support
#' @return Object of class \code{\linkS4class{TrapezoidalFuzzyNumber}}
#' @export
#'
#' @family TrapezoidalFuzzyNumber-method
TrapezoidalFuzzyNumber <- function(a1, a2, a3, a4)
{
   new("TrapezoidalFuzzyNumber", a1=a1, a2=a2, a3=a3, a4=a4)
}


#' @title
#' Creates a Triangular Fuzzy Number
#'
#' @description
#' For convenience, objects of class \code{\linkS4class{TrapezoidalFuzzyNumber}}
#' may be created with this function.
#'
#' @details
#' Currently there is no separate class of a Triangular Fuzzy Number.
#'
#' @param a1 a number specyfing left bound of the support
#' @param amid a number specyfing the core
#' @param a4 a number specyfing right bound of the support
#' @return Object of class \code{\linkS4class{TrapezoidalFuzzyNumber}}
#' @export
#'
#' @family TrapezoidalFuzzyNumber-method
TriangularFuzzyNumber <- function(a1, amid, a4)
{
   new("TrapezoidalFuzzyNumber", a1=a1, a2=amid, a3=amid, a4=a4)
}
