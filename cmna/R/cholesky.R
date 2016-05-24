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

#' @title Cholesky Decomposition
#'
#' @description
#' Decompose a matrix into the Cholesky
#'
#' @param m a matrix
#'
#' @details
#' \code{choleskymatrix} decomposes the matrix \code{m} into the LU
#' decomposition, such that m == L %*% L*.
#'
#' @return the matrix L
#'
#' @family linear
#'
#' @examples
#' (A <- matrix(c(5, 1, 2, 1, 9, 3, 2, 3, 7), 3))
#' (L <- choleskymatrix(A))
#' t(L) %*% L
#'
#' @export
choleskymatrix <- function(m) {
    count.rows <- nrow(m)
    count.cols <- ncol(m)


    L = diag(0, count.rows)
    for(i in 1:count.rows) {
        for(k in 1:i) {
            p.sum <- 0
            for(j in 1:k)
                p.sum <- p.sum + L[j, i] * L[j, k]
            if(i == k)
                L[k, i] <- sqrt(m[i, i] - p.sum)
            else
                L[k, i] <- (m[k, i] - p.sum) / L[k, k]
        }
    }
    return(L)
}
