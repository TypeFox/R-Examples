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

#' @name mcint
#' @rdname mcint
#'
#' @title Monte Carlo Integration
#'
#' @description
#' Simple Monte Carlo Integraton
#'
#' @param f function to integrate
#' @param a the lower-bound of integration
#' @param b the upper-bound of integration
#' @param xdom the domain on \code{x} of integration in two dimensions
#' @param ydom the domain on \code{y} of integration in two dimensions
#' @param m the number of subintervals to calculate
#' @param max.y the largest expected value of the range
#' @param max.z the largest expected value of the range in two dimensions
#'
#' @details
#' The \code{mcint} function uses a simple Monte Carlo algorithm to
#' estimate the value of an integral.  The parameter \code{n} sets the
#' total number of evaluation points.  The parameter \code{max.y} is the
#' maximum expected value of the range of function \code{f}.  The
#' \code{mcint2} provides Monte Carlo integration in two dimensions.
#'
#' @return the value of the integral
#'
#' @family integration
#'
#' @examples
#' f <- function(x) { sin(x)^2 + log(x)}
#' mcint(f, 0, 1)
#' mcint(f, 0, 1, m = 10e6)
#'
#' @importFrom stats runif
#'
#' @export
mcint <- function(f, a, b, m = 1000, max.y = 1) {
    x <- runif(m, min = a, max = b)
    y <- runif(m, min = 0, max = max.y)

    y.hat <- f(x)
    area <- (b - a) * mean(y <= y.hat) * max.y
    return(area)
}

#' @rdname mcint
#' @export
mcint2 <- function(f, xdom, ydom, m = 1000, max.z = 1) {
    xmin <- min(xdom)
    xmax <- max(xdom)
    ymin <- min(ydom)
    ymax <- max(ydom)

    x <- runif(m, min = xmin, max = xmax)
    y <- runif(m, min = ymin, max = ymax)
    z <- runif(m, min = 0, max = max.z)

    z.hat <- f(x, y)
    area <- (xmax - xmin) * (ymax - ymin) *
        mean(z <= z.hat) * max.z
    return(area)
}
