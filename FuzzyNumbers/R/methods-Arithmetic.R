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
#' Arithmetic Operations on Fuzzy Numbers
#'
#' @description
#' Applies arithmetic operations using the extension principle
#' and interval-based calculations.
#'
#' @details
#' Implemented operators: \code{+}, \code{-}, \code{*}, \code{/}
#' for piecewise linear fuzzy numbers.
#' Also some versions may be applied on numeric values and
#' trapezoidal fuzzy numbers.
#'
#' Note that according to the theory the class of PLFNs is not closed
#' under the operations * and /.
#' However, if you operate on a large number of knots,
#' the results should be satisfactory.
#'
#' @param e1 a fuzzy number or single numeric value
#' @param e2 a fuzzy number or single numeric value
#' @return Returns a fuzzy number of the class \linkS4class{PiecewiseLinearFuzzyNumber}
#' or \linkS4class{TrapezoidalFuzzyNumber}.
#'
#' @details
#' Thanks to Jan Caha for suggestions on PLFN operations.
#'
#' @usage
#' \S4method{+}{numeric,FuzzyNumber}(e1, e2) # e2 + e1
#'
#' \S4method{+}{TrapezoidalFuzzyNumber,TrapezoidalFuzzyNumber}(e1, e2)
#'
#' \S4method{+}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#'
#' \S4method{+}{PiecewiseLinearFuzzyNumber,numeric}(e1, e2)
#'
#' \S4method{+}{PiecewiseLinearFuzzyNumber,FuzzyNumber}(e1, e2) # calls as.PiecewiseLinearFuzzyNumber()
#'
#' \S4method{-}{numeric,FuzzyNumber}(e1, e2) # e2*(-1) + e1
#'
#' \S4method{-}{TrapezoidalFuzzyNumber,TrapezoidalFuzzyNumber}(e1, e2)
#'
#' \S4method{-}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#'
#' \S4method{-}{PiecewiseLinearFuzzyNumber,numeric}(e1, e2)
#'
#' \S4method{-}{PiecewiseLinearFuzzyNumber,FuzzyNumber}(e1, e2) # calls as.PiecewiseLinearFuzzyNumber()
#'
#' \S4method{-}{FuzzyNumber,ANY}(e1, e2) # -e1
#'
#' \S4method{*}{numeric,FuzzyNumber}(e1, e2) # e2 * e1
#'
#' \S4method{*}{TrapezoidalFuzzyNumber,numeric}(e1, e2)
#'
#' \S4method{*}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#'
#' \S4method{*}{PiecewiseLinearFuzzyNumber,FuzzyNumber}(e1, e2) # calls as.PiecewiseLinearFuzzyNumber()
#'
#' \S4method{*}{PiecewiseLinearFuzzyNumber,numeric}(e1, e2)
#'
#' \S4method{/}{PiecewiseLinearFuzzyNumber,numeric}(e1, e2)
#'
#' \S4method{/}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#'
#' \S4method{/}{PiecewiseLinearFuzzyNumber,FuzzyNumber}(e1, e2) # calls as.PiecewiseLinearFuzzyNumber()
#'
#' @name Arithmetic
#' @rdname Arithmetic-methods
#' @docType methods
#' @family FuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family TrapezoidalFuzzyNumber-method
#' @family extension_principle
#' @exportMethod /
#' @exportMethod -
#' @exportMethod *
#' @exportMethod +
#' @aliases *,numeric,FuzzyNumber-method
#'          +,numeric,FuzzyNumber-method
#'          -,numeric,FuzzyNumber-method
#'          *,TrapezoidalFuzzyNumber,numeric-method
#'          +,TrapezoidalFuzzyNumber,TrapezoidalFuzzyNumber-method
#'          -,TrapezoidalFuzzyNumber,TrapezoidalFuzzyNumber-method
#'          -,FuzzyNumber,ANY-method
#'          +,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#'          -,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#'          *,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#'          /,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#'          +,PiecewiseLinearFuzzyNumber,FuzzyNumber-method
#'          -,PiecewiseLinearFuzzyNumber,FuzzyNumber-method
#'          /,PiecewiseLinearFuzzyNumber,FuzzyNumber-method
#'          *,PiecewiseLinearFuzzyNumber,FuzzyNumber-method
#'          +,PiecewiseLinearFuzzyNumber,numeric-method
#'          -,PiecewiseLinearFuzzyNumber,numeric-method
#'          *,PiecewiseLinearFuzzyNumber,numeric-method
#'          /,PiecewiseLinearFuzzyNumber,numeric-method
invisible(NULL)





setMethod(
   "*",
   signature(e1 = "numeric", e2 = "FuzzyNumber"),
   function (e1, e2)
   {
      e2*e1
   }
)





setMethod(
   "+",
   signature(e1 = "numeric", e2 = "FuzzyNumber"),
   function (e1, e2)
   {
      e2+e1
   }
)




setMethod(
   "-",
   signature(e1 = "numeric", e2 = "FuzzyNumber"),
   function (e1, e2)
   {
      e2*(-1) + e1
   }
)




setMethod(
   "*",
   signature(e1 = "TrapezoidalFuzzyNumber", e2 = "numeric"),
   function (e1, e2)
   {
      stopifnot(length(e2) == 1);
      TrapezoidalFuzzyNumber(min(e2*e1@a1, e2*e1@a4),
                             min(e2*e1@a2, e2*e1@a3),
                             max(e2*e1@a2, e2*e1@a3),
                             max(e2*e1@a1, e2*e1@a4)
      )
   }
)





setMethod(
   "+",
   signature(e1 = "TrapezoidalFuzzyNumber", e2 = "TrapezoidalFuzzyNumber"),
   function (e1, e2)
   {
      TrapezoidalFuzzyNumber(e1@a1+e2@a1, e1@a2+e2@a2, e1@a3+e2@a3, e1@a4+e2@a4)
   }
)




setMethod(
   "-",
   signature(e1 = "TrapezoidalFuzzyNumber", e2 = "TrapezoidalFuzzyNumber"),
   function (e1, e2)
   {
      TrapezoidalFuzzyNumber(e1@a1-e2@a4, e1@a2-e2@a3, e1@a3-e2@a2, e1@a4-e2@a1)
   }
)






setMethod(
   "-",
   signature(e1 = "FuzzyNumber"),
   function (e1, e2)   # unary minus
   {
      e1*(-1)
   }
)




setMethod(
   "+",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
   function (e1, e2)
   {
      knot.alpha <- unique(sort(c(e1@knot.alpha, e2@knot.alpha)))
      knot.n <- length(knot.alpha)

      if (!isTRUE(all.equal(e1@knot.alpha, knot.alpha)))
         e1 <- as.PiecewiseLinearFuzzyNumber(e1, knot.n=knot.n, knot.alpha=knot.alpha)

      if (!isTRUE(all.equal(e2@knot.alpha, knot.alpha)))
         e2 <- as.PiecewiseLinearFuzzyNumber(e2, knot.n=knot.n, knot.alpha=knot.alpha)

      # using the extension principle and interval-based arithmetic operations
      PiecewiseLinearFuzzyNumber(e1@a1+e2@a1, e1@a2+e2@a2, e1@a3+e2@a3, e1@a4+e2@a4,
                                 knot.n=knot.n, knot.alpha=knot.alpha,
                                 knot.left=e1@knot.left+e2@knot.left,
                                 knot.right=e1@knot.right+e2@knot.right)
   }
)




setMethod(
   "-",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
   function (e1, e2)
   {
      knot.alpha <- unique(sort(c(e1@knot.alpha, e2@knot.alpha)))
      knot.n <- length(knot.alpha)

      if (!isTRUE(all.equal(e1@knot.alpha, knot.alpha)))
         e1 <- as.PiecewiseLinearFuzzyNumber(e1, knot.n=knot.n, knot.alpha=knot.alpha)

      if (!isTRUE(all.equal(e2@knot.alpha, knot.alpha)))
         e2 <- as.PiecewiseLinearFuzzyNumber(e2, knot.n=knot.n, knot.alpha=knot.alpha)

      # using the extension principle and interval-based arithmetic operations
      PiecewiseLinearFuzzyNumber(knot.alpha=knot.alpha,
                                 knot.left=c(e1@a1,e1@knot.left,e1@a2)-rev(c(e2@a3,e2@knot.right,e2@a4)),
                                 knot.right=c(e1@a3,e1@knot.right,e1@a4)-rev(c(e2@a1,e2@knot.left,e2@a2))
      )
   }
)






setMethod(
   "*",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
   function (e1, e2)
   {
      knot.alpha <- unique(sort(c(e1@knot.alpha, e2@knot.alpha)))
      knot.n <- length(knot.alpha)

      if (!isTRUE(all.equal(e1@knot.alpha, knot.alpha)))
         e1 <- as.PiecewiseLinearFuzzyNumber(e1, knot.n=knot.n, knot.alpha=knot.alpha)

      if (!isTRUE(all.equal(e2@knot.alpha, knot.alpha)))
         e2 <- as.PiecewiseLinearFuzzyNumber(e2, knot.n=knot.n, knot.alpha=knot.alpha)

      e1l <- c(e1@a1,e1@knot.left,e1@a2)
      e2l <- c(e2@a1,e2@knot.left,e2@a2)
      e1r <- rev(c(e1@a3,e1@knot.right,e1@a4))
      e2r <- rev(c(e2@a3,e2@knot.right,e2@a4))

      p1 <- e1l*e2l
      p2 <- e1l*e2r
      p3 <- e1r*e2l
      p4 <- e1r*e2r
      
      # warning("Piecewise linear fuzzy numbers are not closed under the `*' operation.")

      # using the extension principle and interval-based arithmetic operations
      PiecewiseLinearFuzzyNumber(knot.alpha=knot.alpha,
                                 knot.left=pmin(p1, p2, p3, p4),
                                 knot.right=rev(pmax(p1, p2, p3, p4))
      )
   }
)




setMethod(
   "/",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
   function (e1, e2)
   {
      if (e2@a1 <= 0 && e2@a4 >= 0)
         stop("division by zero")

      knot.alpha <- unique(sort(c(e1@knot.alpha, e2@knot.alpha)))
      knot.n <- length(knot.alpha)

      if (!isTRUE(all.equal(e1@knot.alpha, knot.alpha)))
         e1 <- as.PiecewiseLinearFuzzyNumber(e1, knot.n=knot.n, knot.alpha=knot.alpha)

      if (!isTRUE(all.equal(e2@knot.alpha, knot.alpha)))
         e2 <- as.PiecewiseLinearFuzzyNumber(e2, knot.n=knot.n, knot.alpha=knot.alpha)

      e1l <- c(e1@a1,e1@knot.left,e1@a2)
      e2r <- 1.0/c(e2@a1,e2@knot.left,e2@a2)
      e1r <- rev(c(e1@a3,e1@knot.right,e1@a4))
      e2l <- 1.0/rev(c(e2@a3,e2@knot.right,e2@a4))

      p1 <- e1l*e2l
      p2 <- e1l*e2r
      p3 <- e1r*e2l
      p4 <- e1r*e2r

      # warning("Piecewise linear fuzzy numbers are not closed under the `/' operation.")

      # using the extension principle and interval-based arithmetic operations
      PiecewiseLinearFuzzyNumber(knot.alpha=knot.alpha,
                                 knot.left=pmin(p1, p2, p3, p4),
                                 knot.right=rev(pmax(p1, p2, p3, p4))
      )
   }
)






setMethod(
   "+",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "FuzzyNumber"),
   function (e1, e2)
   {
      e1 + as.PiecewiseLinearFuzzyNumber(e2, knot.n=e1@knot.n, knot.alpha=e1@knot.alpha)
   })



setMethod(
   "-",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "FuzzyNumber"),
   function (e1, e2)
   {
      e1 - as.PiecewiseLinearFuzzyNumber(e2, knot.n=e1@knot.n, knot.alpha=e1@knot.alpha)
   })



setMethod(
   "*",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "FuzzyNumber"),
   function (e1, e2)
   {
      e1 * as.PiecewiseLinearFuzzyNumber(e2, knot.n=e1@knot.n, knot.alpha=e1@knot.alpha)
   })



setMethod(
   "/",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "FuzzyNumber"),
   function (e1, e2)
   {
      e1 / as.PiecewiseLinearFuzzyNumber(e2, knot.n=e1@knot.n, knot.alpha=e1@knot.alpha)
   })




setMethod(
   "+",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "numeric"),
   function (e1, e2)
   {
      e1 + as.PiecewiseLinearFuzzyNumber(e2, knot.n=e1@knot.n, knot.alpha=e1@knot.alpha)
   })



setMethod(
   "-",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "numeric"),
   function (e1, e2)
   {
      e1 - as.PiecewiseLinearFuzzyNumber(e2, knot.n=e1@knot.n, knot.alpha=e1@knot.alpha)
   })



setMethod(
   "*",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "numeric"),
   function (e1, e2)
   {
      stopifnot(length(e2) == 1)
      kl <-     c(e1@a1, e1@knot.left,  e1@a2)
      kr <- rev(c(e1@a3, e1@knot.right, e1@a4))
      kmin <- pmin(e2*kl, e2*kr)
      kmax <- pmax(e2*kl, e2*kr)

      PiecewiseLinearFuzzyNumber(
         kmin[1],
         kmin[length(kmin)],
         kmax[length(kmax)],
         kmax[1],
         knot.n=e1@knot.n,
         knot.alpha=e1@knot.alpha,
         knot.left= kmin[-c(1,length(kmin))],
         knot.right=rev(kmax[-c(1,length(kmax))])
      )
   }
)



setMethod(
   "/",
   signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "numeric"),
   function (e1, e2)
   {
      e1 / as.PiecewiseLinearFuzzyNumber(e2, knot.n=e1@knot.n, knot.alpha=e1@knot.alpha)
   })
