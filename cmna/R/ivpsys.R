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

#' @title Initial value problems for systems of ordinary differential equations
#'
#' @name ivpsys
#' @rdname ivpsys
#'
#' @description
#' solve initial value problems for systems ordinary differential equations
#'
#' @param f function to integrate
#' @param x0 the initial value of x
#' @param y0 the vector initial values of y
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

#' @rdname ivpsys
#' @export
eulersys <- function(f, x0, y0, h, n) {
    x <- x0
    y <- y0

    ## If y0 values are named, the data frame names them!
    ## The value names produced by f(x, y) should match.
    values <- data.frame(x = x, t(y0))
    for(i in 1:n) {
        y0 <- y0 + h * f(x0, y0)
        x0 <- x0 + h
        values <- rbind(values, data.frame(x = x0, t(y0)))
    }

    return(values)
}
