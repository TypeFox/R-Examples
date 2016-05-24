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
#' Get Basic Information on a Fuzzy Number in a String
#'
#' @description
#' This method is especially useful if you would like to
#' generate LaTeX equations defining a fuzzy numbers.
#'
#' @details
#' Consider calling the \code{\link{cat}} function on the resulting
#' string.
#'
#' Thanks to Jan Caha for suggesting the \code{toLaTeX} arg.
#'
#' @param x a fuzzy number
#' @param toLaTeX logical; should LaTeX code be output?
#' @param varnameLaTeX character; variable name to be included in equations
#' 
#' @return Returns a character vector.
#'
#' @usage
#' \S4method{as.character}{FuzzyNumber}(x, toLaTeX=FALSE, varnameLaTeX="A")
#'
#' \S4method{as.character}{TrapezoidalFuzzyNumber}(x, toLaTeX=FALSE, varnameLaTeX="A")
#'
#' \S4method{as.character}{PiecewiseLinearFuzzyNumber}(x, toLaTeX=FALSE, varnameLaTeX="A")
#'
#' \S4method{as.character}{PowerFuzzyNumber}(x, toLaTeX=FALSE, varnameLaTeX="A")
#'
#'
#' @exportMethod as.character
#' @name as.character
#' @docType methods
#' @rdname as.character-methods
#' @family FuzzyNumber-method
#' @family TriangularFuzzyNumber-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @family PowerFuzzyNumber-method
#' @family conversion
#' @aliases as.character,FuzzyNumber-method
#'          as.character,PiecewiseLinearFuzzyNumber-method
#'          as.character,TrapezoidalFuzzyNumber-method
#'          as.character,PowerFuzzyNumber-method
invisible(NULL)


setMethod(
   f="as.character",
   signature(x="FuzzyNumber"),
   definition=function(x, toLaTeX=FALSE, varnameLaTeX="A") {
      if (identical(toLaTeX, FALSE)) {
         sprintf("Fuzzy number with:\n   support=[%g,%g],\n      core=[%g,%g].\n",
                 x@a1, x@a4, x@a2, x@a3)
      }
      else {
         res <- sprintf(paste(
            "\\[",
            "\\mu_{%s}(x) = \\left\\{",
            "\\begin{array}{lll}",
            "0      & \\text{for} & x\\in(-\\infty,%g), \\\\",
            "l_{%s}(x) & \\text{for} & x\\in[%g,%g), \\\\",
            "1      & \\text{for} & x\\in[%g,%g], \\\\",
            "r_{%s}(x) & \\text{for} & x\\in(%g,%g], \\\\",
            "0      & \\text{for} & x\\in(%g,+\\infty), \\\\",
            "\\end{array}",
            "\\right.",
            "\\]",
            sep="\n"),
            varnameLaTeX,
            x@a1,
            varnameLaTeX, x@a1, x@a2,
            x@a2, x@a3,
            varnameLaTeX, x@a3, x@a4,
            x@a4
         )

         res <- paste(res, sprintf(paste(
            "where $l_{%s}=\\mathtt{left}_A((x%+g)/%g)$,",
            "$r_{%s}=\\mathtt{right}_A((x%+g)/%g)$.",
            sep="\n"),
           varnameLaTeX, -x@a1, x@a2-x@a1,
           varnameLaTeX, -x@a3, x@a4-x@a3
         ), sep="\n")

         res <- paste(res, "\n", sep="")

         res <- paste(res, sprintf(paste(
            "\\[",
            "{%s}_\\alpha = [{%s}_L(\\alpha), {%s}_U(\\alpha)],",
            "\\]",
            "where ${%s}_L(\\alpha)=%g%+g\\,\\mathtt{lower}_{%s}(\\alpha)$,",
            "${%s}_U(\\alpha)=%g%+g\\,\\mathtt{upper}_{%s}(\\alpha)$.",
            sep="\n"),
           varnameLaTeX, varnameLaTeX, varnameLaTeX,
           varnameLaTeX, x@a1, x@a2-x@a1, varnameLaTeX,
           varnameLaTeX, x@a3, x@a4-x@a3, varnameLaTeX
         ), sep="\n")

         res
      }
   }
)



setMethod(
   f="as.character",
   signature(x="PiecewiseLinearFuzzyNumber"),
   definition=function(x, toLaTeX=FALSE, varnameLaTeX="A") {
      if (identical(toLaTeX, FALSE)) {
         sprintf("Piecewise linear fuzzy number with %g knot(s),\n   support=[%g,%g],\n      core=[%g,%g].\n",
            x@knot.n, x@a1, x@a4, x@a2, x@a3)
      }
      else {

         a <- c(0, x@knot.alpha, 1)
         l <- c(x@a1, x@knot.left, x@a2)
         r <- c(x@a3, x@knot.right, x@a4)
         n <- x@knot.n

         res <- sprintf(paste(
            "\\[",
            "\\mu_{%s}(x) = \\left\\{",
            "\\begin{array}{lll}",
            "0      & \\text{for} & x\\in(-\\infty,%g), \\\\",
            sep="\n"),
            varnameLaTeX,
            x@a1)



         for (i in 1:(n+1))
            res <- paste(res, sprintf("%g%+g\\,(x%+g)/%g & \\text{for} & x\\in[%g,%g), \\\\",
                                      a[i], a[i+1]-a[i], l[i], l[i+1]-l[i], l[i], l[i+1]), sep="\n")

         res <- paste(res,
                      sprintf("1      & \\text{for} & x\\in[%g,%g], \\\\", x@a2, x@a3),
                      sep="\n")

         for (i in 1:(n+1))
            res <- paste(res, sprintf("%g%+g\\,(%g-x)/%g & \\text{for} & x\\in[%g,%g), \\\\",
                                      a[n-i+2], a[n-i+3]-a[n-i+2], r[i+1], r[i+1]-r[i], r[i], r[i+1]), sep="\n")

         res <- paste(res, sprintf(paste(
            "0      & \\text{for} & x\\in(%g,+\\infty). \\\\",
            "\\end{array}",
            "\\right.",
            "\\]",
            sep="\n"),
            x@a4), sep="\n")

         res <- paste(res, "\n", sep="")

         res <- paste(res, sprintf(paste(
            "\\[",
            "{%s}_\\alpha = [{%s}_L(\\alpha), {%s}_U(\\alpha)],",
            "\\]",
            "where",
            sep="\n"),
            varnameLaTeX, varnameLaTeX, varnameLaTeX
         ), sep="\n")

         res <- paste(res, sprintf(paste(
            "\\[",
            "{%s}_L(\\alpha) = \\left\\{",
            "\\begin{array}{lll}",
            sep="\n"),
            varnameLaTeX,
            x@a1
         ), sep="\n")

         for (i in 1:(n+1))
            res <- paste(res, sprintf("%g%+g\\,(\\alpha%+g)/%g & \\text{for} & \\alpha\\in[%g,%g], \\\\",
                                      l[i], l[i+1]-l[i], -a[i], a[i+1]-a[i], a[i], a[i+1]), sep="\n")

         res <- paste(res,
            "\\end{array}",
            "\\right.",
            "\\]",
            sep="\n")

         res <- paste(res, sprintf(paste(
            "\\[",
            "{%s}_U(\\alpha) = \\left\\{",
            "\\begin{array}{lll}",
            sep="\n"),
              varnameLaTeX,
              x@a1
         ), sep="\n")

         for (i in 1:(n+1))
            res <- paste(res, sprintf("%g%+g\\,(%g-\\alpha)/%g & \\text{for} & \\alpha\\in[%g,%g]%s \\\\",
                                      r[n-i+2], r[n-i+3]-r[n-i+2], a[i+1], a[i+1]-a[i], a[i], a[i+1],
                                    if(i < n+1) "," else "."), sep="\n")

         res <- paste(res,
             "\\end{array}",
             "\\right.",
             "\\]",
          sep="\n")

         res
      }
   }
)



setMethod(
   f="as.character",
   signature(x="TrapezoidalFuzzyNumber"),
   definition=function(x, toLaTeX=FALSE, varnameLaTeX="A") {
      if (identical(toLaTeX, FALSE)) {
         sprintf("Trapezoidal fuzzy number with:\n   support=[%g,%g],\n      core=[%g,%g].\n",
                 x@a1, x@a4, x@a2, x@a3)
      }
      else {
         res <- sprintf(paste(
            "\\[",
            "\\mu_{%s}(x) = \\left\\{",
            "\\begin{array}{lll}",
            "0      & \\text{for} & x\\in(-\\infty,%g), \\\\",
            "(x%+g)/%g & \\text{for} & x\\in[%g,%g), \\\\",
            "1      & \\text{for} & x\\in[%g,%g], \\\\",
            "(%g-x)/%g & \\text{for} & x\\in(%g,%g], \\\\",
            "0      & \\text{for} & x\\in(%g,+\\infty). \\\\",
            "\\end{array}",
            "\\right.",
            "\\]",
            sep="\n"),
            varnameLaTeX,
            x@a1,
            -x@a1, x@a2-x@a1, x@a1, x@a2,
            x@a2, x@a3,
            x@a4, x@a4-x@a3, x@a3, x@a4,
            x@a4
         )

         res <- paste(res, "\n", sep="")

         res <- paste(res, sprintf(paste(
            "\\[",
            "{%s}_\\alpha = [%g%+g\\,\\alpha, %g%+g\\,\\alpha].",
            "\\]",
            sep="\n"),
            varnameLaTeX,
            x@a1, x@a2-x@a1,
            x@a4, -(x@a4-x@a3)
         ), sep="\n")

         res
      }
   }
)



setMethod(
   f="as.character",
   signature(x="PowerFuzzyNumber"),
   definition=function(x, toLaTeX=FALSE, varnameLaTeX="A") {
      if (identical(toLaTeX, FALSE)) {
         sprintf("Fuzzy number given by power functions, and:\n   support=[%g,%g],\n      core=[%g,%g].\n",
                 x@a1, x@a4, x@a2, x@a3)
      }
      else {
         res <- sprintf(paste(
            "\\[",
            "\\mu_{%s}(x) = \\left\\{",
            "\\begin{array}{lll}",
            "0      & \\text{for} & x\\in(-\\infty,%g), \\\\",
            "((x%+g)/%g)^{%g} & \\text{for} & x\\in[%g,%g), \\\\",
            "1      & \\text{for} & x\\in[%g,%g], \\\\",
            "((%g-x)/%g)^{%g} & \\text{for} & x\\in(%g,%g], \\\\",
            "0      & \\text{for} & x\\in(%g,+\\infty), \\\\",
            "\\end{array}",
            "\\right.",
            "\\]",
            sep="\n"),
            varnameLaTeX,
            x@a1,
            -x@a1, x@a2-x@a1, x@p.left, x@a1, x@a2,
            x@a2, x@a3,
            x@a4, x@a4-x@a3, x@p.right, x@a3, x@a4,
            x@a4
         )

         res <- paste(res, "\n", sep="")

         res <- paste(res, sprintf(paste(
            "\\[",
            "{%s}_\\alpha = [%g%+g\\,\\alpha^{%g}, %g%+g\\,(1-\\alpha^{%g})].",
            "\\]",
            sep="\n"),
                                   varnameLaTeX,
                                   x@a1, x@a2-x@a1, 1/x@p.left,
                                   x@a3, x@a4-x@a3, 1/x@p.right
         ), sep="\n")

         res
      }
   }
)
