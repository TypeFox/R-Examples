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

#' @title Hill climbing
#'
#' @name hillclimbing
#' @rdname hillclimbing
#'
#' @description
#' Use hill climbing to find the global minimum
#'
#' @param f function representing the derivative of \code{f}
#' @param x an initial estimate of the minimum
#' @param h the step size
#' @param m the maximum number of iterations
#'
#' @details
#'
#' Hill climbing
#'
#' @return the \code{x} value of the minimum found
#'
#' @family optimz
#'
#' @examples
#' f <- function(x) {
#'     (x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2
#' }
#' hillclimbing(f, c(0,0))
#' hillclimbing(f, c(-1,-1))
#' hillclimbing(f, c(10,10))
#'
#' @importFrom stats runif
#' @importFrom stats rnorm

#' @export
hillclimbing <- function(f, x, h = 1, m = 1e3) {
    n <- length(x)

    xcurr <- x
    ycurr <- f(x)

    for(i in 1:m) {
        xnext <- xcurr
        i <- ceiling(runif(1, 0, n))
        xnext[i] <- rnorm(1, xcurr[i], h)
        ynext <- f(xnext)
        if(ynext < ycurr) {
            xcurr <- xnext
            ycurr <- ynext
        }
    }

    return(xcurr)
}

