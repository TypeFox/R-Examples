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

#' @title Secant Method
#'
#' @description
#' The secant method for root finding
#'
#' @param f function to integrate
#' @param x an initial estimate of the root
#' @param tol the error tolerance
#' @param m the maximum number of iterations
#'
#' @details
#'
#' The secant method for root finding extends Newton's method to
#' estimate the derivative.  It will return when the interval between
#' them is less than \code{tol}, the error tolerance.  However, this
#' implementation also stop if after \code{m} iterations.
#'
#' @return the real root found
#'
#' @family optimz
#'
#' @examples
#' f <- function(x) { x^3 - 2 * x^2 - 159 * x - 540 }
#' secant(f, 1)
#'
#' @export
secant <- function(f, x, tol = 1e-3, m = 100) {
    i <- 0

    oldx <- x
    oldfx <- f(x)
    x <- oldx + 10 * tol

    while(abs(x - oldx) > tol) {
        i <- i + 1
        if (i > m)
            stop("No solution found")

        fx <- f(x)
        newx <- x - fx * ((x - oldx) / (fx - oldfx))
        oldx <- x
        oldfx <- fx
        x <- newx
    }

    return(x)
}
