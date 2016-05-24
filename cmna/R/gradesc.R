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

#' @title Gradient descent
#'
#' @name gradient
#' @rdname gradient
#'
#' @description
#' Use gradient descent to find local minima
#'
#' @param fp function representing the derivative of \code{f}
#' @param x an initial estimate of the minima
#' @param h the step size
#' @param tol the error tolerance
#' @param m the maximum number of iterations
#'
#' @details
#'
#' Gradient descent can be used to find local minima of functions.  It
#' will return an approximation based on the step size \code{h} and
#' \code{fp}.  The \code{tol} is the error tolerance, \code{x} is the
#' initial guess at the minimum.  This implementation also stops after
#' \code{m} iterations.
#'
#' @return the \code{x} value of the minimum found
#'
#' @family optimz
#'
#' @examples
#' fp <- function(x) { x^3 + 3 * x^2 - 1 }
#' graddsc(fp, 0)
#'
#' f <- function(x) { (x[1] - 1)^2 + (x[2] - 1)^2 }
#' fp <-function(x) {
#'     x1 <- 2 * x[1] - 2
#'     x2 <- 8 * x[2] - 8
#'
#'     return(c(x1, x2))
#' }
#' gd(fp, c(0, 0), 0.05)

#' @export
graddsc <- function(fp, x, h = 1e-3, tol = 1e-4, m = 1e3) {
    iter <- 0

    oldx <- x
    x = x - h * fp(x)

    while(abs(x - oldx) > tol) {
        iter <- iter + 1
        if(iter > m)
            stop("No solution found")
        oldx <- x
        x = x - h * fp(x)
    }

    return(x)
}

#' @rdname gradient
#' @export
gradasc <- function(fp, x, h = 1e-3, tol = 1e-4, m = 1e3) {
    iter <- 0

    oldx <- x
    x = x + h * fp(x)

    while(abs(x - oldx) > tol) {
        iter <- iter + 1
        if(iter > m)
            stop("No solution found")
        oldx <- x
        x = x + h * fp(x)
    }

    return(x)
}

#' @rdname gradient
#' @export
gd <- function(fp, x, h = 1e2, tol = 1e-4, m = 1e3) {
    iter <- 0

    oldx <- x
    x = x - h * fp(x)

    while(vecnorm(x - oldx) > tol) {
        iter <- iter + 1
        if(iter > m)
            return(x)
        oldx <- x
        x = x - h * fp(x)
    }

    return(x)
}
