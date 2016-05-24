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

#' @name horner
#' @rdname horner
#'
#' @title Horner's rule
#'
#' @description
#' Use Horner's rule to evaluate a polynomial
#'
#' @param x a vector of x values to evaluate the polynomial
#' @param betas vector of coefficients of x
#'
#' @details
#' This function implements Horner's rule for fast polynomial
#' evaluation.  The implementation expects \code{x} to be a vector of x
#' values at which to evaluate the polynomial. The parameter \code{betas}
#' is a vector of coefficients of \emph{x}.  The vector order is such
#' that the first element is the constant term, the second element is
#' the coefficient of \emph{x}, the so forth to the highest degreed
#' term.  Terms with a 0 coefficient should have a 0 element in the
#' vector.
#'
#' The function \code{rhorner} implements the the Horner algorithm
#' recursively.
#'
#' The function \code{naivepoly} implements a polynomial evaluator using
#' the straightforward algebraic approach.
#'
#' The function \code{betterpoly} implements a polynomial evaluator using
#' the straightforward algebraic approach with cached \emph{x} terms.
#'
#' @return the value of the function at \code{x}
#'
#' @family algebra
#'
#' @examples
#' b <- c(2, 10, 11)
#' x <- 5
#' horner(x, b)
#' b <- c(-1, 0, 1)
#' x <- c(1, 2, 3, 4)
#' horner(x, b)
#' rhorner(x, b)

#' @export
horner <- function(x, betas) {
    y <- rep(0, length(x))

    for(i in length(betas):1)
        y <- betas[i] + x * y

    return(y)
}

#' @rdname horner
#' @export
rhorner <- function(x, betas) {
    n <- length(betas)

    if(n == 1)
        return(betas)

    return(betas[1] + x * rhorner(x, betas[2:n]))
}

#' @rdname horner
#' @export
naivepoly <- function(x, betas) {
    y <- rep(0, length(x))

    for(i in 1:length(betas))
        y <- y + betas[i] * (x ^ (i - 1))

    return(y)
}

#' @rdname horner
#' @export
betterpoly <- function(x, betas) {
    y <- rep(0, length(x))
    cached.x <- 1

    for(i in 1:length(betas)) {
        y <- y + betas[i] * cached.x
        cached.x <- cached.x * x
    }

    return(y)
}
