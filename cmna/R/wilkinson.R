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

#' @title Wilkinson's Polynomial
#'
#' @description
#' Wilkinson's polynomial
#'
#' @param x the \code{x}-value
#' @param w the number of terms in the polynomial
#'
#' @details
#' Wilkinson's polynomail is a terrible joke played on numerical
#' analysis.  By tradition, the function is f(x) = (x - 1)(x - 2)...(x -
#' 20), giving a function with real roots at each integer from 1 to 20.
#' This function is generalized and allows for \code{n} and the function
#' value is f(x) = (x - 1)(x - 2)...(x - n).  The default of \code{n} is
#' 20.
#'
#' @return the value of the function at \code{x}
#'
#' @family polynomials
#'
#' @examples
#' wilkinson(0)
#'
#' @export
wilkinson <- function(x, w = 20) {

    if(w == 1)
        return(x - 1)
    return((x - w) * wilkinson(x, w - 1))
}
