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

#' @title Boundary value problems
#'
#' @name bvp
#' @rdname bvp
#'
#' @description
#' solve boundary value problems for ordinary differential equations
#'
#' @param x proposed initial \code{x}-value
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
#' bvpexample(-2)
#' bvpexample(-1)
#' bvpexample(0)
#' bvpexample(1)
#' bvpexample(2)
#' ## (bvp.b <- bisection(bvpexample, 0, 1))
#' ## (bvp.s <- secant(bvpexample, 0))
#'
#' @importFrom utils tail
#'
#' @rdname bvp
#' @export
bvpexample <- function(x) {
    x0 <- 0
    y0 <- c(y1 = 1, y2 = x)
    yn <- 1

    odesystem <- function(x, y) {
        y1 <- y[2]
        y2 <- y[1]^2 - 2

        return(c(y1 = y1, y2 = y2))
    }

    z <- eulersys(odesystem, x0, y0, 1 / 1000, 1000)
    tail(z$y1, 1) - yn
}

#' @rdname bvp
#' @export
bvpexample10 <- function(x) {
  x0 <- 0
  y0 <- c(y1 = 1, y2 = x)
  yn <- 1

  odesystem <- function(x, y) {
      y1 <- y[2]
      y2 <- y[1]^2 - 2

      return(c(y1 = y1, y2 = y2))
  }

  z <- eulersys(odesystem, x0, y0, 1/ 10, 10)
  tail(z$y1, 1) - yn
}
