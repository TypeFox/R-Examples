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

#' @title Adaptive Integration
#'
#' @description
#' Adaptive integration
#'
#' @param f function to integrate
#' @param a the a-bound of integration
#' @param b the b-bound of integration
#' @param n the maximum recursive depth
#' @param tol the maximum error tolerance
#'
#' @details
#' The \code{adaptint} function uses Romberg's rule to calculate the
#' integral of the function \code{f} over the interval from \code{a}
#' to \code{b}.  The parameter \code{n} sets the number of intervals
#' to use when evaluating.  Additional options are passed to the
#' function \code{f} when evaluating.
#'
#' @return the value of the integral
#'
#' @family integration
#' @family newton-cotes
#' @family adaptive
#'
#' @examples
#' f <- function(x) { sin(x)^2 + log(x) }
#' adaptint(f, 1, 10, n = 4)
#' adaptint(f, 1, 10, n = 5)
#' adaptint(f, 1, 10, n = 10)
#'
#' @export
adaptint <- function(f, a, b, n = 10, tol = 1e-6) {
    if(n == 1)
        area <- midpt(f, a, b, m = 2)
    else {
        q1 <- midpt(f, a, b, m = 1)
        q2 <- midpt(f, a, b, m = 2)
        if(abs(q1 - q2) > 3 * tol) {
            n = n - 1
            tol <- tol / 2
            c <- (a + b) / 2
            lt <- adaptint(f, a, c, n = n, tol = tol)
            rt <- adaptint(f, c, b, n = n, tol = tol)
            area <- lt + rt
        } else
            area <- q2
    }

    return(area)
}
