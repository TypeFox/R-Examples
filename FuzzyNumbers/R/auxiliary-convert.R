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
#' Convert a Given Side Function to Side Generating Function
#'
#' @description
#' The resulting function linearly scales the input
#' and passes it to the original function.
#'
#' @details
#' The function works for x1<x2 and x1>x2.
#'
#'
#' @param f a function defined on [x1,x2]
#' @param x1 numeric vector of length 1; if longer, only the first element is used
#' @param x2 numeric vector of length 1; if longer, only the first element is used
#' @return Returns a new function defined on [0,1] (scaled input).
#'
#' @seealso \code{\linkS4class{FuzzyNumber}}
#' @family auxiliary
#' @export
convertSide <- function(f, x1, x2)
{
   stopifnot(is.function(f) && is.numeric(x1) && is.numeric(x2))
   e <- new.env()
   e$f <- f
   e$x1 <- x1[1]
   e$dx <- dx <- x2[1]-x1[1] # dx <- ... to avoid CHECK problems
   stopifnot(e$dx != 0 && is.finite(e$x1) && is.finite(e$dx))
   side <- function(x) {
      f(x1+dx*x)
   }
   environment(side) <- e
   side
}


#' @title
#' Convert a Given Upper/Lower Alpha-Cut Function to an Alpha-Cut Generating Function
#'
#' @description
#' The resulting function calls the original function and then
#' linearly scales its output.
#'
#' @param f a function into [y1,y2]
#' @param y1 numeric vector of length 1
#' @param y2 numeric vector of length 1
#' @return Returns a new function defined on [0,1] (scaled input).
#'
#' @seealso \code{\linkS4class{FuzzyNumber}}
#' @family auxiliary
#' @export
convertAlpha <- function(f, y1, y2)
{
   stopifnot(is.function(f) && is.numeric(y1) && is.numeric(y2))
   e <- new.env()
   e$f <- f
   e$y1 <- y1[1]
   e$dy <- dy <- y2[1]-y1[1] # dy <- ... to avoid CHECK problems
   stopifnot(e$dy != 0 && is.finite(e$y1) && is.finite(e$dy))
   alpha <- function(x) {
      (f(x)-y1)/dy
   }
   environment(alpha) <- e
   alpha
}
