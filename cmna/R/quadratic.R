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

#' @rdname quadratic
#' @name quadratic
#'
#' @title The quadratic equation.
#'
#' @description
#' Find the zeros of a quadratic equation.
#'
#' @param b2 the coefficient of the x^2 term
#' @param b1 the coefficient of the x term
#' @param b0 the constant term
#'
#' @details
#' \code{quadratic} and \code{quadratic2} implement the quadratic
#' equation from standard algebra in two different ways.  The
#' \code{quadratic} function is susceptible to cascading numerical error
#' and the \code{quadratic2} has reduced potential error.
#'
#' @return numeric vector of solutions to the equation
#'
#' @family algebra
#'
#' @examples
#' quadratic(1, 0, -1)
#' quadratic(4, -4, 1)
#' quadratic2(1, 0, -1)
#' quadratic2(4, -4, 1)

#' @export
quadratic <- function(b2, b1, b0) {
    t1 <- sqrt(b1^2 - 4 * b2 * b0)
    t2 <- 2 * b2

    x1 <- - (b1 + t1) / t2
    x2 <- - (b1 - t1) / t2
    return(c(x1, x2))
}

#' @rdname quadratic
#' @export
quadratic2 <- function(b2, b1, b0) {
    t1 <- sqrt(b1^2 - 4 * b2 * b0)
    t2 <- 2 * b0

    x1 <- t2 / (-b1 - t1)
    x2 <- t2 / (-b1 + t1)

    ## Reverse the order so they come
    ## back the same as quadratic()
    return(c(x2, x1))
}
