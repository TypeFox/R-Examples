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

#' @title Invert a matrix
#'
#' @description
#' Invert the matrix using Gaussian elimination
#'
#' @param m a matrix
#'
#' @details
#' \code{invmatrix} invertse the given matrix using Gaussian elimination
#' and returns the result.
#'
#' @return the inverted matrix
#'
#' @family linear
#'
#' @examples
#' A <- matrix(c(1, 2, -7, -1, -1, 1, 2, 1, 5), 3)
#' refmatrix(A)
#'
#' @export
invmatrix <- function(m) {
    count.rows <- nrow(m)
    count.cols <- ncol(m)

    if(count.rows != count.cols)
        stop("Matrix must be square")

    I = diag(count.cols)
    tmp <- cbind(m, I)
    tmp <- rrefmatrix(tmp)
    res <- tmp[, (count.cols + 1):(count.cols * 2)]
    return(res)
}
