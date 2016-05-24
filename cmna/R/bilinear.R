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

#' @title Bilinear interpolation
#'
#' @description
#' Finds a bilinear interpolation bounded by four points
#'
#' @param x vector of two x values representing \code{x_1} and \code{x_2}
#' @param y vector of two y values representing \code{y_1} and \code{y_2}
#' @param z 2x2 matrix if \code{z} values
#' @param newx vector of new \code{x} values to interpolate
#' @param newy vector of new \code{y} values to interpolate
#'
#' @details
#' \code{bilinear} finds a bilinear interpolation bounded by four corners
#'
#' @return a vector of interpolated z values at (\code{x}, \code{y})
#'
#' @family interp
#' @family algebra
#'
#' @examples
#' x <- c(2, 4)
#' y <- c(4, 7)
#' z <- matrix(c(81, 84, 85, 89), nrow = 2)
#' newx <- c(2.5, 3, 3.5)
#' newy <- c(5, 5.5, 6)
#' bilinear(x, y, z, newx, newy)
#'
#' @export
bilinear <- function(x, y, z, newx, newy) {

    ## Find intermediate values along the x-axis, first
    z1 <- (x[2] - newx) * z[1,1] + (newx - x[1]) * z[1,2]
    z1 <- z1 / (x[2] - x[1])
    z2 <- (x[2] - newx) * z[2,1] + (newx - x[1]) * z[2,2]
    z2 <- z2 / (x[2] - x[1])

    ## Then interpolate along the y-axis
    z <- (y[2] - newy) * z1 + (newy - y[1]) * z2
    z <- z / (y[2] - y[1])

    return(z)
}

