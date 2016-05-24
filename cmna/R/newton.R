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

#' @title Newton's method
#'
#' @description
#' Use Newton's method to find real roots
#'
#' @param f function to integrate
#' @param fp function representing the derivative of \code{f}
#' @param x an initial estimate of the root
#' @param tol the error tolerance
#' @param m the maximum number of iterations
#'
#' @details
#'
#' Newton's method finds real roots of a function, but requires knowing
#' the function derivative.  It will return when the interval between
#' them is less than \code{tol}, the error tolerance.  However, this
#' implementation also stops after \code{m} iterations.
#'
#' @return the real root found
#'
#' @family optimz
#'
#' @examples
#' f <- function(x) { x^3 - 2 * x^2 - 159 * x - 540 }
#' fp <- function(x) {3 * x^2 - 4 * x - 159 }
#' newton(f, fp, 1)
#'
#' @export
newton <- function(f, fp, x, tol = 1e-3, m = 100) {
    iter <- 0

    oldx <- x
    x <- oldx + 10 * tol

    while(abs(x - oldx) > tol) {
        iter <- iter + 1
        if(iter > m)
            stop("No solution found")
        oldx <- x
        x <- x - f(x) / fp(x)
    }

    return(x)
}
