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

#' @title Initial value problems
#'
#' @name ivp
#' @rdname ivp
#'
#' @description
#' solve initial value problems for ordinary differential equations
#'
#' @param f function to integrate
#' @param x0 the initial value of x
#' @param y0 the initial value of y
#' @param h selected step size
#' @param n the number of steps
#'
#' @details
#' The \code{euler} method implements the Euler method for solving
#' differential equations.  The code{midptivp} method solves initial
#' value problems using the second-order Runge-Kutta method.  The
#' \code{rungekutta4} method is the fourth-order Runge-Kutta method.
#'
#' @return a data frame of \code{x} and \code{y} values
#'
#' @examples
#' f <- function(x, y) { y / (2 * x + 1) }
#' ivp.euler <- euler(f, 0, 1, 1/100, 100)
#' ivp.midpt <- midptivp(f, 0, 1, 1/100, 100)
#' ivp.rk4 <- rungekutta4(f, 0, 1, 1/100, 100)

#' @export
euler <- function(f, x0, y0, h, n) {
    x <- x0
    y <- y0

    for(i in 1:n) {
        y0 <- y0 + h * f(x0, y0)
        x0 <- x0 + h
        x <- c(x, x0)
        y <- c(y, y0)
    }

    return(data.frame(x = x, y = y))
}

#' @rdname ivp
#' @export
midptivp <- function(f, x0, y0, h, n) {
    x <- x0
    y <- y0

    for(i in 1:n) {
        s1 <- h * f(x0, y0)
        s2 <- h * f(x0 + h / 2, y0 + s1 / 2)
        y0 <- y0 + s2

        x0 <- x0 + h
        x <- c(x, x0)
        y <- c(y, y0)
    }

    return(data.frame(x = x, y = y))
}

#' @rdname ivp
#' @export
rungekutta4 <- function(f, x0, y0, h, n) {
    x <- x0
    y <- y0

    for(i in 1:n) {
        s1 <- h * f(x0, y0)
        s2 <- h * f(x0 + h / 2, y0 + s1 / 2)
        s3 <- h * f(x0 + h / 2, y0 + s2 / 2)
        s4 <- h * f(x0 + h, y0 + s3)
        y0 <- y0 + s1 / 6 + s2 / 3 + s3 / 3 + s4 / 6

        x0 <- x0 + h
        x <- c(x, x0)
        y <- c(y, y0)
    }

    return(data.frame(x = x, y = y))
}

#' @rdname ivp
#' @export
adamsbashforth <- function(f, x0, y0, h, n) {

    ## Quick Euler the value of x1, y1
    y1 <- y0 + h * f(x0, y0)
    x1 <- x0 + h

    x <- c(x0, x1)
    y <- c(y0, y1)
    n <- n - 1

    for(i in 1:n) {
        yn <- y1 + 1.5 * h * f(x1, y1) - .5 * h * f(x0, y0)
        xn <- x1 + h

        y0 <- y1
        x0 <- x1
        y1 <- yn
        x1 <- xn

        y <- c(y, y1)
        x <- c(x, x1)
    }

    return(data.frame(x = x, y = y))
}
