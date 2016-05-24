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

#' @title Wave Equation using
#'
#' @name wave
#' @rdname wave
#'
#' @description
#' solve heat equation via forward-time central-space method
#'
#' @param u the initial values of u
#' @param alpha the thermal diffusivity coefficient
#' @param xdelta the change in \code{x} at each step in \code{u}
#' @param tdelta the time step
#' @param n the number of steps to take
#'
#' @details
#' The \code{heat} solves the heat equation using the forward-time
#' central-space method in one-dimension.
#'
#' @return a matrix of u values at each time step
#'
#' @examples
#' speed <- 2
#' x0 <- 0
#' xdelta <- .05
#' x <- seq(x0, 1, xdelta)
#' m <- length(x)
#' u <- sin(x * pi * 2)
#' u[11:21] <- 0
#' tdelta <- .02
#' n <- 40
#' z <- wave(u, speed, xdelta, tdelta, n)

#' @export
wave <- function(u, alpha, xdelta, tdelta, n) {
    m <- length(u)
    uarray <- matrix(u, nrow = 1)
    newu <- u

    h <- ((alpha * tdelta) / (xdelta))^2

    ## Initial the zeroth timestep
    oldu <- rep(0, m)
    oldu[2:(m - 1)] <- u[2:(m - 1)] + h *
        (u[1:(m - 2)] - 2 * u[2:(m - 1)] + u[3:m]) / 2

    ## Now iterate
    for(i in 1:n) {
        ustep1 <- (2 * u - oldu)
        ustep2 <- u[1:(m - 2)] - 2 * u[2:(m - 1)] + u[3:m]
        newu <- ustep1 + h * c(0, ustep2, 0)
        oldu <- u
        u <- newu
        uarray <- rbind(uarray, u)
    }

    return(uarray)
}
