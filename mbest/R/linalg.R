# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


proj.psd <- function(x)
{
    modified <- FALSE
    e <- eigen(x, symmetric=TRUE)
    l <- e$values
    u <- e$vectors

    if (any(l < 0)) {
        xproj <- u %*% diag(ifelse(l > 0, l, 0), nrow(x)) %*% t(u)
        dimnames(xproj) <- dimnames(x)
        x <- xproj
        modified <- TRUE
    }

    attr(x, "modified") <- modified
    x
}


pseudo.solve <- function(a, b)
{
    withCallingHandlers({
        a.chol <- chol(a, pivot=TRUE)
    }, warning = function(w) {
        if (conditionMessage(w) == "the matrix is either rank-deficient or indefinite")
            invokeRestart("muffleWarning")
    })

    rank <- attr(a.chol, "rank")
    n <- nrow(a.chol)

    pivot <- attr(a.chol, "pivot")[seq_len(rank)]
    R <- a.chol[seq_len(rank), seq_len(rank)]

    deficient <- rank < n

    if (!deficient) {
        if (missing(b)) {
            x <- matrix(0, n, n)
            x[pivot,pivot] <- chol2inv(R)
        } else {
            if (is.matrix(b)) {
                y <- b[pivot,,drop=FALSE]
                x <- matrix(0, n, ncol(b))
                x[pivot,] <- backsolve(R, backsolve(R, y, transpose=TRUE))
            } else  {
                y <- b[pivot]
                x <- numeric(n)
                x[pivot] <- backsolve(R, backsolve(R, y, transpose=TRUE))
            }
        }
    } else {
        a.eigen <- eigen(a, symmetric=TRUE)
        u <- a.eigen$vectors[,seq_len(rank),drop=FALSE]
        l <- a.eigen$values[seq_len(rank)]
        if (missing(b)) {
            x <- u %*% (t(u) / l)
        } else {
            x <- u %*% ((t(u) %*% b) / l)
        }
    }

    attr(x, "deficient") <- deficient
    x
}

