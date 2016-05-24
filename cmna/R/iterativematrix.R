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

#' @rdname iterativematrix
#' @name iterativematrix
#'
#' @title Solve a matrix using iterative methods
#'
#' @description
#' Solve a matrix using iterative methods.
#'
#' @param A a square matrix representing the coefficients of a linear
#' system
#' @param b a vector representing the right-hand side of the linear
#' system
#' @param tol is a number representing the error tolerence
#' @param maxiter is the maximum number of iterations
#'
#' @details
#' \code{jacobi} finds the solution using Jacobi iteration.
#' Jacobi iteration depends on the matrix being diagonally-dominate.
#' The tolerence is specified the norm of the solution vector.
#'
#' \code{gaussseidel} finds the solution using Gauss-Seidel iteration.
#' Gauss-Seidel iteration depends on the matrix being either
#' diagonally-dominate or symmetric and positive definite.
#'
#' \code{cgmmatrix} finds the solution using the conjugate gradient
#' method.  The conjugate gradient method depends on the matrix being
#' symmetric and positive definite.
#'
#' @return the solution vector
#'
#' @family linear
#'
#' @examples
#' A <- matrix(c(5, 2, 1, 2, 7, 3, 3, 4, 8), 3)
#' b <- c(40, 39, 55)
#' jacobi(A, b)
#'
#' @export
jacobi <- function(A, b, tol = 10e-7, maxiter = 100) {
    n <- length(b)
    iter <- 0

    Dinv <- diag(1 / diag(A))
    R <- A - diag(diag(A))
    x <- rep(0, n)
    newx <- rep(tol, n)

    while(vecnorm(newx - x) > tol & iter < maxiter) {
        x <- newx
        newx <- Dinv %*% (b - R %*% x)
        iter <- iter + 1
    }

    return(as.vector(x))
}

#' @rdname iterativematrix
#' @export
gaussseidel <- function(A, b, tol = 10e-7, maxiter = 100) {
    n <- length(b)
    iter <- 0

    L <- U <- A
    L[upper.tri(A, diag = FALSE)] <- 0
    U[lower.tri(A, diag = TRUE)] <- 0
    Linv <- solve(L)

    x <- rep(0, n)
    newx <- rep(tol * 10, n)

    while(vecnorm(newx - x) > tol & iter < maxiter) {
        x <- newx
        newx <- Linv %*% (b - U %*% x)
        iter <- iter + 1
    }

    return(as.vector(x))
}

#' @rdname iterativematrix
#' @export
cgmmatrix <- function(A, b, tol = 10e-7, maxiter = 100) {
    n <- length(b)
    iter <- 0

    x <- rep(0, n)
    newx <- rep(tol * 10, n)

    p <- r <- b - A %*% x
    while(vecnorm(r) > tol & iter < maxiter) {
        a <- as.numeric((t(r) %*% r) / t(p) %*% A %*% p)
        newx <- x + a * p
        newr <- r - a * A %*% p
        beta <- as.numeric(t(newr) %*% newr / (t(r) %*% r))
        p <- newr + beta * p
        r <- newr
        x <- newx
        iter <- iter + 1
    }

    return(as.vector(x))
}
