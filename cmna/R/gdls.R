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

#' @name gdls
#'
#' @title Least squares with graident descent
#'
#' @description
#' Solve least squares with graident descent
#'
#' @param A a square matrix representing the coefficients of a linear system
#' @param b a vector representing the right-hand side of the linear system
#' @param alpha the learning rate
#' @param tol the expected error tolerance
#' @param m the maximum number of iterations
#'
#' @details
#'
#' \code{gdls} solves a linear system using gradient descent.
#'
#' @return the modified matrix
#'
#' @family linear
#'
#' @examples
#' head(b <- iris$Sepal.Length)
#' head(A <- matrix(cbind(1, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width), ncol = 4))
#' gdls(A, b, alpha = 0.05, m = 10000)
#'
#' @export
gdls <- function(A, b, alpha = 0.05, tol = 1e-6, m = 1e5) {
    iter <- 0
    n <- ncol(A)
    theta <- matrix(rep(0, n))
    oldtheta = theta + 10 * tol

    while(vecnorm(oldtheta - theta) > tol) {
        if((iter <- iter + 1) > m)
            return(theta)
        e <- (A %*% theta - b)
        d <- (t(A) %*% e) / length(b)
        oldtheta <- theta
        theta <- theta - alpha * d
    }

    return(theta)
}
