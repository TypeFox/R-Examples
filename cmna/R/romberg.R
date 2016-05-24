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

#' @title Romberg Integration
#'
#' @description
#' Romberg's adaptive integration
#'
#' @param f function to integrate
#' @param a the lowerbound of integration
#' @param b the upperbound of integration
#' @param m the maximum number of iterations
#' @param tab if \code{TRUE}, return the table of values
#'
#' @details
#' The \code{romberg} function uses Romberg's rule to calculate the
#' integral of the function \code{f} over the interval from \code{a}
#' to \code{b}.  The parameter \code{m} sets the number of intervals
#' to use when evaluating.  Additional options are passed to the
#' function \code{f} when evaluating.
#'
#' @return the value of the integral
#'
#' @family integration
#' @family newton-cotes
#'
#' @examples
#' f <- function(x) { sin(x)^2 + log(x)}
#' romberg(f, 1, 10, m = 3)
#' romberg(f, 1, 10, m = 5)
#' romberg(f, 1, 10, m = 10)
#'
#' @export
romberg <- function(f, a, b, m, tab = FALSE) {
    R <- matrix(NA, nrow = m, ncol = m)

    R[1, 1] <- trap(f, a, b, m = 1)
    for(j in 2:m) {
        R[j, 1] <- trap(f, a, b, m = 2^(j - 1))
        for(k in 2:j) {
            k4 <- 4^(k - 1)
            R[j, k] <- k4 * R[j, k - 1] - R[j - 1, k - 1]
            R[j, k] <- R[j, k] / (k4 - 1)
        }
    }

    if(tab == TRUE)
        return(R)
    return(R[m, m])
}
