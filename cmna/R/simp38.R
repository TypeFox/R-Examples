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

#' @title Simpson's 3/8 rule
#'
#' @description
#' Use Simpson's 3/8 rule to integrate a function
#'
#' @param f function to integrate
#' @param a the a-bound of integration
#' @param b the b-bound of integration
#' @param m the number of subintervals to calculate
#'
#' @details
#' The \code{simp38} function uses Simpson's 3/8 rule to calculate the
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
#' f <- function(x) { sin(x)^2 + log(x) }
#' simp38(f, 1, 10, m = 10)
#' simp38(f, 1, 10, m = 100)
#' simp38(f, 1, 10, m = 1000)
#'
#' @export
simp38 <- function(f, a, b, m = 100) {
    x.ends = seq(a, b, length.out = m + 1)
    y.ends = f(x.ends)
    x.midh = (2 * x.ends[2:(m + 1)] + x.ends[1:m]) / 3
    x.midl = (x.ends[2:(m + 1)] + 2 * x.ends[1:m]) / 3
    y.midh = f(x.midh)
    y.midl = f(x.midl)

    p.area = sum(y.ends[2:(m+1)] + 3 * y.midh[1:m] + 3 *
                     y.midl[1:m] + y.ends[1:m])
    p.area = p.area * abs(b - a) / (8 * m)
    return(p.area)
}
