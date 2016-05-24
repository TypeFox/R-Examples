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

#' @title LU Decomposition
#'
#' @description
#' Decompose a matrix into lower- and upper-triangular matrices
#'
#' @param m a matrix
#'
#' @details
#' \code{lumatrix} decomposes the matrix \code{m} into the LU
#' decomposition, such that m == L %*% U.
#'
#' @return list with matrices L and U representing the LU decomposition
#'
#' @family linear
#'
#' @examples
#' A <- matrix(c(1, 2, -7, -1, -1, 1, 2, 1, 5), 3)
#' lumatrix(A)
#'
#' @export
lumatrix <- function(m) {
    count.rows <- nrow(m)
    count.cols <- ncol(m)
    piv <- 1

    P <- L <- diag(count.cols)
    for(row.curr in 1:count.rows) {
        if(piv <= count.cols) {
            i <- row.curr
            while(m[i, piv] == 0 && i < count.rows) {
                i <- i + 1
                if(i > count.rows) {
                    i <- row.curr
                    piv <- piv + 1
                    if(piv > count.cols)
                        return(list(P = P, L = L, U = m))
                }
            }
            if(i != row.curr) {
                m <- swaprows(m, i, row.curr)
                P <- swaprows(P, i, row.curr)
            }
            for(j in row.curr:count.rows)
                if(j != row.curr) {
                    k <- m[j, piv] / m[row.curr, piv]
                    m <- replacerow(m, row.curr, j, -k)
                    L[j, piv] <- k
                }
            piv <- piv + 1
        }
    }

    return(list(P = P, L = L, U = m))
}
