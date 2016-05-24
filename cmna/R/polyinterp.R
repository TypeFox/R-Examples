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

#' @title Polynomial interpolation
#'
#' @description
#' Finds a polynomial function interpolating the given points
#'
#' @param x a vector of x values
#' @param y a vector of y values
#'
#' @details
#' \code{polyinterp} finds a polynomial that interpolates the given points.
#'
#' @return a polynomial equation's coefficients
#'
#' @family interp
#' @family algebra
#'
#' @examples
#' x <- c(1, 2, 3)
#' y <- x^2 + 5 * x - 3
#' f <- polyinterp(x, y)
#'
#' @export
polyinterp <- function(x, y) {
    if(length(x) != length(y))
        stop("Length of x and y vectors must be the same")

    n <- length(x) - 1
    vandermonde <- rep(1, length(x))
    for(i in 1:n) {
        xi <- x^i
        vandermonde <- cbind(xi, vandermonde)
    }
    beta <- solve(vandermonde, y)

    names(beta) <- NULL
    return(rev(beta))
}
