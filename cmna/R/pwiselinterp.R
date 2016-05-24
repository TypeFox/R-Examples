## Copyright (c) 2016, James P. Howard, II <jh@jameshoward.us>
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' @title Piecewise linear interpolation
#'
#' @description
#' Finds a piecewise linear function that interpolates the data points
#'
#' @param x a vector of x values
#' @param y a vector of y values
#'
#' @details
#' \code{pwiselinterp} finds a piecewise linear function that
#' interpolates the data points.  For each x-y ordered pair, there
#' function finds the unique line interpolating them.  The function will
#' return a data.frame with three columns.
#'
#' The column \code{x} is the upper bound of the domain for the given
#' piece.  The columns \code{m} and \code{b} represent the coefficients
#' from the y-intercept form of the linear equation, y = mx + b.
#'
#' The matrix will contain length(x) rows with the first row having m
#' and b of NA.
#'
#' @return a matrix with the linear function components
#'
#' @family interp
#' @family algebra
#'
#' @examples
#' x <- c(5, 0, 3)
#' y <- c(4, 0, 3)
#' f <- pwiselinterp(x, y)
#'
#' @export
pwiselinterp <- function(x, y) {
    n <- length(x) - 1

    y <- y[order(x)]
    x <- x[order(x)]

    mvec <- bvec <- c()

    for(i in 1:n) {
        p <- linterp(x[i], y[i], x[i + 1], y[i + 1])
        mvec <- c(mvec, p[2])
        bvec <- c(bvec, p[1])
    }

    return(list(m = mvec, b = bvec))
}
