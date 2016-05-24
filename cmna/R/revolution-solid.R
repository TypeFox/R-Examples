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

#' @name revolution-solid
#' @rdname revolution-solid
#'
#' @title Volumes of solids of revolution
#'
#' @description
#' Find the volume of a solid of revolution
#'
#' @param f function of revolution
#' @param a lower-bound of the solid
#' @param b upper-bound of the solid
#'
#' @details
#'
#' The functions \code{discmethod} and \code{shellmethod} implement the
#' algorithms for finding the volume of solids of revolution.  The
#' \code{discmethod} function is suitable for volumes revolved around
#' the \code{x}-axis and the \code{shellmethod} function is suitable for
#' volumes revolved around the \code{y}-axis.
#'
#' @return the volume of the solid
#'
#' @family integration
#'
#' @examples
#' f <- function(x) { x^2 }
#' shellmethod(f, 1, 2)
#' discmethod(f, 1, 2)
#'
#' @export
shellmethod <- function(f, a, b) {
    solid <- function(x) { return(x * f(x)) }

    return(2 * pi * trap(solid, a, b))
}

#' @rdname revolution-solid
#' @export
discmethod <- function(f, a, b) {
    solid <- function(x) { return(pi * (f(x))^2) }

    return(midpt(solid, a, b))
}
