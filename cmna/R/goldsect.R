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

#' @title Golden Section Search
#'
#' @name goldsect
#' @rdname goldsect
#'
#' @description
#' Use golden section search to find local extrema
#'
#' @param f function to integrate
#' @param a the a bound of the search region
#' @param b the b bound of the search region
#' @param tol the error tolerance
#' @param m the maximum number of iterations
#'
#' @details
#'
#' The golden section search  method functions by repeatedly dividing the interval
#' between \code{a} and \code{b} and will return when the
#' interval between them is less than \code{tol}, the error tolerance.
#' However, this implementation also stop if after \code{m}
#' iterations.
#'
#' @return the \code{x} value of the minimum found
#'
#' @family optimz
#'
#' @examples
#' f <- function(x) { x^2 - 3 * x + 3 }
#' goldsectmin(f, 0, 5)
#'
#' @export
goldsectmin <- function(f, a, b, tol = 1e-3, m = 100) {
    iter <- 0
    phi <- (sqrt(5) - 1) / 2

    a.star <- b - phi * abs(b - a)
    b.star <- a + phi * abs(b - a)

    while (abs(b - a) > tol) {
        iter <- iter + 1
        if (iter > m)
            break

        if(f(a.star) < f(b.star)) {
            b <- b.star
            b.star <- a.star
            a.star <- b - phi * abs(b - a)
        } else {
            a <- a.star
            a.star <- b.star
            b.star <- a + phi * abs(b - a)
        }
    }

    return((a + b) / 2)
}

#' @rdname goldsect
#' @export
goldsectmax <- function(f, a, b, tol = 1e-3, m = 100) {
    iter <- 0
    phi <- (sqrt(5) - 1) / 2

    a.star <- b - phi * abs(b - a)
    b.star <- a + phi * abs(b - a)

    while (abs(b - a) > tol) {
        iter <- iter + 1
        if (iter > m)
            break

        if(f(a.star) > f(b.star)) {
            b <- b.star
            b.star <- a.star
            a.star <- b - phi * abs(b - a)
        } else {
            a <- a.star
            a.star <- b.star
            b.star <- a + phi * abs(b - a)
        }
    }

    return((a + b) / 2)
}
