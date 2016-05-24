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

#' @title Simulated annealing
#'
#' @name sa
#' @rdname sa
#'
#' @description
#' Use simulated annealing to find the global minimum
#'
#' @param f function representing \code{f}
#' @param x an initial estimate of the minimum
#' @param temp the initial temperature
#' @param rate the cooling rate
#'
#' @details
#'
#' Simulated annealing finds a global minimum by mimicing the
#' metallurgical process of annealing.
#'
#' @return the \code{x} value of the minimum found
#'
#' @family optimz
#'
#' @examples
#' f <- function(x) { x^6 - 4 * x^5 - 7 * x^4 + 22 * x^3 + 24 * x^2 + 2}
#' sa(f, 0)
#'
#' f <- function(x) { (x[1] - 1)^2 + (x[2] - 1)^2 }
#' sa(f, c(0, 0), 0.05)
#'
#' @importFrom stats runif
#' @importFrom stats rnorm
#'
#' @rdname sa
#' @export
sa <- function(f, x, temp = 1e4, rate = 1e-4) {
    step = 1 - rate
    n <- length(x)

    xbest <- xcurr <- xnext <- x
    ybest <- ycurr <- ynext <- f(x)

    while(temp > 1) {
        temp <- temp * step
        i <- ceiling(runif(1, 0, n))
        xnext[i] <- rnorm(1, xcurr[i], temp)
        ynext <- f(xnext)
        accept <- exp(-(ynext - ycurr) / temp)
        if(ynext < ycurr || runif(1) < accept) {
            xcurr <- xnext
            ycurr <- ynext
        }
        if(ynext < ybest) {
            xbest <- xcurr
            ybest <- ycurr
        }
    }

    return(xbest)
}

#' @rdname sa
#' @export
tspsa <- function(x, temp = 1e2, rate = 1e-4) {
    step = 1 - rate
    n <- nrow(x)

    xbest <- xcurr <- xnext <- c(1:n)
    ynext <- 0
    for(i in 2:n) {
        a <- xnext[i - 1]
        b <- xnext[i]
        ynext <- ynext + vecnorm(x[a,] - x[b,])
    }
    a <- xnext[1]
    b <- xnext[n]
    ynext <- ynext + vecnorm(x[a,] - x[b,])
    ybest <- ycurr <- ynext

    while(temp > 1) {
        temp <- temp * step
        i <- ceiling(runif(1, 1, n))
        xnext <- xcurr
        temporary <- xnext[i]
        xnext[i] <- xnext[i - 1]
        xnext[i - 1] <- temporary
        ynext <- 0
        for(i in 2:n) {
          a <- xnext[i - 1]
          b <- xnext[i]
          ynext <- ynext + vecnorm(x[a,] - x[b,])
        }
        a <- xnext[1]
        b <- xnext[n]
        ynext <- ynext + vecnorm(x[a,] - x[b,])
        accept <- exp(-(ynext - ycurr) / temp)
        if(ynext < ycurr || runif(1) < accept) {
            xcurr <- xnext
            ycurr <- ynext
        }
        if(ynext < ybest) {
            xbest <- xcurr
            ybest <- ycurr
        }
    }
    return(list(order = xbest, distance = ybest))
}
