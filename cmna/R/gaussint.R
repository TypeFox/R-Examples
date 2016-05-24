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

#' @title Gaussian integration method driver
#'
#' @description
#' Use the Gaussian method to evaluate integrals
#'
#' @param f function to integrate
#' @param m number of evaluation points
#' @param x list of evaluation points
#' @param w list of weights
#'
#' @details
#' The \code{gaussint} function uses the Gaussian integration to
#' evaluate an integral.  The function itself is a driver and expects
#' the integration points and associated weights as options.
#'
#' @return the value of the integral
#'
#' @family integration
#'
#' @examples
#' w = c(1, 1)
#' x = c(-1 / sqrt(3), 1 / sqrt(3))
#' f <- function(x) { x^3 + x + 1 }
#' gaussint(f, x, w)
#'
#' @export
gaussint <- function(f, x, w) {
    y <- f(x)

    return(sum(y * w))
}

#' @rdname gaussint
#' @export
gauss.legendre <- function(f, m = 5) {
    p <- paste("gauss.legendre.", m, sep = "")
    params <- eval(parse(text = p))

    return(gaussint(f, params$x, params$w))
}

#' @rdname gaussint
#' @export
gauss.laguerre <- function(f, m = 5) {
  p <- paste("gauss.laguerre.", m, sep = "")
  params <- eval(parse(text = p))

  return(gaussint(f, params$x, params$w))
}

#' @rdname gaussint
#' @export
gauss.hermite <- function(f, m = 5) {
  p <- paste("gauss.hermite.", m, sep = "")
  params <- eval(parse(text = p))

  return(gaussint(f, params$x, params$w))
}
