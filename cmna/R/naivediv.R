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

#' @name division
#' @rdname division
#'
#' @title Algorithms for divisions
#'
#' @description
#' Algorithms for division that provide a quotient and remainder.
#'
#' @param m the dividend
#' @param n the divisor
#'
#' @details
#' The \code{naivediv} divides \code{m} by \code{n} by using repeated
#' division.  The \code{longdiv} function uses the long division
#' algorithm in binary.
#'
#' @return the quotient and remainder as a list
#'
#' @family algebra
#'
#' @examples
#' a <- floor(runif(1, 1, 1000))
#' b <- floor(runif(1, 1, 100))
#' naivediv(a, b)
#' longdiv(a, b)
#'
#' @export
naivediv <- function(m, n) {
    quot <- 0
    r <- m

    if(n == 0)
        stop("Attempted division by 0")

    while(r >= n) {
        quot <- quot + 1
        r <- r - n
    }

    return(list(quot = quot, r = r))
}

#' @rdname division
#' @export
longdiv <- function(m, n) {
    quot <- 0
    r <- 0

    if(n == 0)
        stop("Attempted division by 0")

    for(i in 31:0) {
        r <- bitwShiftL(r, 1)
        r <- r + bitwAnd(bitwShiftR(m, i), 1)
        if(r >= n) {
            r <- r - n
            quot <- quot + bitwShiftL(1, i)
        }
    }

    return(list(quot = quot, r = r))
}
