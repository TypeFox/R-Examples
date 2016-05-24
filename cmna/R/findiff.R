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

#' @name findiff
#' @rdname findiff
#'
#' @title Finite Differences
#'
#' @description
#' Finite differences formulas
#'
#' @param f function to differentiate
#' @param x the \code{x}-value to differentiate at
#' @param h the step-size for evaluation
#' @param tol the error tolerance for \code{symdiff}
#' @param m the maximum number of convergence steps in \code{symdiff}
#' @param n the maximum number of convergence steps in \code{rdiff}
#'
#' @details
#'
#' The \code{findiff} formula uses the finite differences formula to
#' find the derivative of \code{f} at \code{x}.  The value of \code{h}
#' is the step size of the evaluation. The function \code{findiff2}
#' provides the second derivative.
#'
#' @return the value of the derivative
#'
#' @family differentiation
#'
#' @examples
#' findiff(sin, pi, 1e-3)
#' symdiff(sin, pi, 1e-3)
#'
#' @export
findiff <- function(f, x, h) {
    return((f(x + h) - f(x)) / h)
}

#' @rdname findiff
#' @export
symdiff <- function(f, x, h = tol * 10,
                    tol = 1e-3, m = 100) {
    i <- 0

    lastdx <- (f(x + h) - f(x - h)) / (2 * h)
    while(i < m) {
        h <- h / 2
        dx <- (f(x + h) - f(x - h)) / (2 * h)
        if((abs(dx - lastdx) < tol))
            return(dx)
        lastdx <- dx
    }
}

#' @rdname findiff
#' @export
findiff2 <- function(f, x, h) {
    return((f(x + h) - 2 * f(x) + f(x - h)) / h^2)
}

#' @rdname findiff
#' @export
rdiff <- function(f, x, n = 10, h = 1e-4) {
    if(n == 1)
        return(symdiff(f, x, h = h))

    dx <- (4 * rdiff(f, x, n = n - 1, h = h / 2) -
               symdiff(f, x, h = h)) / 3
    return(dx)
}
