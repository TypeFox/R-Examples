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

#' @title Nearest interpolation
#'
#' @description
#' Find the nearest neighbor for a set of data points
#'
#' @param p matrix of variable values, each row is a data point
#' @param y vector of values, each entry corresponds to one row in \code{p}
#' @param q vector of variable values, each entry corresponds to one column of \code{p}
#'
#' @details
#' \code{nn} finds the n-dimensional nearest neighbor for given datapoint
#'
#' @return an interpolated value for \var{q}
#'
#' @family interp
#'
#' @examples
#' p <- matrix(floor(runif(100, 0, 9)), 20)
#' y <- floor(runif(20, 0, 9))
#' q <- matrix(floor(runif(5, 0, 9)), 1)
#' nn(p, y, q)
#'
#' @export
nn <- function(p, y, q) {
    if(ncol(p) != ncol(q))
        stop("p and q must have same number of columns")

    ## Repeat the rows of q to simplfy the  calculation
    qprime <- t(matrix(rep(q, nrow(p)), ncol(p)))
    d <- sqrt(rowSums((p - qprime)^2))
    return(y[which.min(d)])
}
