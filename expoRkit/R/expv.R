### R interface for calling Expokit Fortran routines for computing the
### matrix exponential. This file contains S4-methods for different matrix
### classes from the Matrix package and SparseM package together with a
### method for ordinary matrices. 
###
###     Copyright (C) 2012 Niels Richard Hansen.
###
### This program is free software; you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation; either version 2, or (at your option) any
### later version.
###
### These functions are distributed in the hope that they will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, a copy is available at
### http://www.r-project.org/Licenses/

##' Methods for computing the solution of the ODE \deqn{w'(t) = x w(t)
##' + u} with initial condition \eqn{w(0) = v} at one or more time
##' points. 
##'
##' Analytically the solution is given as \deqn{w(t) = \exp(tx)v + t
##' \phi(tx)u}{w(t) = exp(tx)v + tphi(tx)u} with \eqn{\phi(z) =
##' (\exp(z)-1)/z}{phi(z) = (exp(z)-1)/z}. For large matrices \eqn{x}
##' the computation of the full matrices \eqn{\exp(tx)}{exp(tx)} and
##' \eqn{\phi(tx)}{phi(tx)} is slow. An alternative is to compute the
##' solution \eqn{w} directly for a given initial value \eqn{v}
##' using Krylov subspace methods. This is, in particular, efficient
##' for large sparse matrices.
##'
##' Note that if \eqn{Q} is a rate matrix for a homogeneous continuous
##' time Markov process (non-negative off-diagonals and row-sums 0)
##' and \eqn{v} is a probability vector (the initial distribution at
##' time point 0), then the distribution at time point \eqn{t} solves
##' \deqn{w'(t) = Q^T w(t).} In this case we want to take \eqn{x} to
##' be \eqn{Q^T} to get the desired solution.
##'
##' The solution is computed using the Fortran package
##' Expokit. Methods are available for matrix classes implemented in
##' the Matrix package as well as the SparseM package. The
##' implementation avoids the computation of the full matrix
##' exponential of \eqn{tx} and the approach is advantageous when we
##' want to compute \eqn{w(t)} for one or a few initial values
##' \eqn{v}. The full matrix exponential should \emph{not} be computed
##' this way by looping over \eqn{n} different initial values.
##' 
##' Though there is a method implemented for ordinary (dense)
##' matrices, such a matrix is simply coerced into a
##' \code{CsparseMatrix} before the solution is computed. It is
##' recommended that large sparse matrices are stored and handled as
##' such, e.g. using the classes and methods from the Matrix
##' package. Dense intermediates should be avoided.
##'
##' The \code{x} matrix is allowed to be a dense complex matrix in
##' which case \code{v} and \code{u} are also allowed to be complex.
##' @title Matrix exponentiation and more. 
##' @param x a matrix. 
##' @param v a \code{numeric} vector. The initial value. 
##' @param t a \code{numeric} vector of time points at which the
##' solution is computed.
##' @param u a \code{numeric} vector. Default \code{NULL}.
##' @param Markov \code{logical}. If \code{TRUE} the matrix is taken
##' to be an rate matrix and steps are taken to ensure that the
##' computed result is a probability vector. Default \code{FALSE}.
##' @param transpose \code{logical}. If \code{TRUE} transpose the
##' matrix before the solution is computed. Default equals
##' \code{Markov}.
##' @param ... other arguments passed to \code{\link{Rexpv}}.
##' @return An \eqn{n} by \eqn{k} matrix with \eqn{k} the length of
##' \code{t} and \eqn{n} the dimension of the matrix \code{x}. The \eqn{i}'th
##' column contains the solution of the ODE at time point \eqn{i}.
##' @author Niels Richard Hansen \email{Niels.R.Hansen@@math.ku.dk}
##' @seealso \code{\link{Rexpv}}, \code{\link[Matrix]{expm}},
##' \code{\link[expm]{expm}}
##' @examples
##' ### Small 4 by 4 example.
##' x <- matrix(c(-1, 0, 1, 0,
##'               0, -2, 0, -3,
##'               1, 0, 0, 0,
##'               0, 2, -1, 0),
##'             4, 4)
##' v <- c(1, 1, 1, 1)
##' 
##' require(Matrix)
##' require(SparseM)
##' 
##' w <- cbind(padm(x) %*% v,
##'            expv(x, v),
##'            expv(Matrix(x, sparse = TRUE), v),
##'            expv(as.matrix.coo(x), v),
##'            expv(as.matrix.csr(x), v),
##'            expv(as.matrix.csc(x), v)
##'            )
##' 
##' stopifnot(all.equal(w[, 1], w[, 2]),
##'           all.equal(w[, 1], w[, 3]),
##'           all.equal(w[, 1], w[, 4]),
##'           all.equal(w[, 1], w[, 5]),
##'           all.equal(w[, 1], w[, 6]))
##' 
##' u <- c(2, 0, 1, 1)
##' ex <- padm(x)
##' w <- cbind(ex %*% v + (ex - diag(1, 4)) %*% solve(x, u),
##'            expv(x, v, u = u),
##'            expv(Matrix(x, sparse = TRUE), v, u = u),
##'            expv(as.matrix.coo(x), v, u = u),
##'            expv(as.matrix.csr(x), v, u = u),
##'            expv(as.matrix.csc(x), v, u = u)
##'            )
##' 
##' stopifnot(all.equal(w[, 1], w[, 2]),
##'           all.equal(w[, 1], w[, 3]),
##'           all.equal(w[, 1], w[, 4]),
##'           all.equal(w[, 1], w[, 5]),
##'           all.equal(w[, 1], w[, 6]))
##' 
##' ############################################################
##' ### Linear birth-death Markov process with immigration
##' ############################################################
##'
##' alpha <- 2.1  ## Death rate per individual
##' beta <- 2     ## Birth rate per individual
##' delta <- 20   ## Immigration rate
##' 
##' n <- 500L     ## state space {0, ..., n-1}
##' i <- seq(1, n)
##' rates <- c(alpha * i[-1], ## subdiagonal
##'            -(c(0, alpha * i[-1]) +
##'              c(delta, beta * i[-c(1,n)] + delta, 0)), ## diagonal
##'            c(delta, beta * i[-c(1, n)] + delta)) ## superdiagonal
##' j <- c(i[-n], i, i[-1])
##' i <- c(i[-1], i, i[-n])
##'
##' ## Sparse rate matrix constructed without dense intermediate
##' require(Matrix)
##' Q <- sparseMatrix(i = i, j = j, x = rates, dims = c(n, n))
##'
##' ## Evolution of uniform initial distribution
##' p0 <- rep(1, n)/n
##' time <- seq(0, 10, 0.2)
##' Pt <- expv(Q, p0, t = time, Markov = TRUE)
##' 
##' \dontrun{
##' matplot(1:n, Pt, type = "l")
##' image(time, 0:(n-1), -t(Pt), col = terrain.colors(100))
##' }
##' @export
##' @docType methods
##' @rdname expv-methods
##' @exportMethod expv
setGeneric("expv", function(x, v, t = 1.0, u = NULL,
                            Markov = FALSE, transpose = Markov, ...)
           standardGeneric("expv")
           )

##' @rdname expv-methods
##' @aliases expv,matrix,vector-method
##' @importFrom Matrix Matrix
setMethod("expv", signature("matrix", "vector"),
          function(x, v, t = 1.0, u = NULL,
                   Markov = FALSE, transpose = Markov, ...) {
            if (is.numeric(x)) {
            ## For a numeric 'matrix' argument coerce the argument into a
            ## sparse matrix in CCS format.
            w <- callGeneric(Matrix(x, sparse = TRUE), v = v, t = t, u = u,
                             Markov = Markov, transpose = transpose, ...)
          } else if (is.complex(x)) {
            ## For a complex 'matrix' the coercion to a sparse matrix
            ## is done manually as the Matrix package does not support
            ## complex matrices (yet).
            pattern <- x != 0
            sparsePattern <- as(pattern, 'Matrix')
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x[pattern],
                              ia = sparsePattern@i + 1,
                              ja = sparsePattern@p + 1,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = FALSE,
                              storage = 'CCS',
                              ...)
              v <- w[, k]
            }
          } else {
            stop("Argument 'x' must be a 'numeric' or 'complex' matrix.")
          }
            w
          }
          )

##' @rdname expv-methods
##' @aliases expv,CsparseMatrix,vector-method
##' @importClassesFrom Matrix CsparseMatrix
##' @importMethodsFrom Matrix t
setMethod("expv", signature("CsparseMatrix", "vector"),
          function(x, v, t = 1.0, u = NULL,
                   Markov = FALSE, transpose = Markov, ...) {
            callGeneric(as(x, "dgCMatrix"), v = v, t = t, u = u,
                        Markov = Markov, transpose = transpose, ...)
        }
        )

##' @rdname expv-methods
##' @aliases expv,dgCMatrix,vector-method
##' @importClassesFrom Matrix dgCMatrix
setMethod("expv", signature("dgCMatrix", "vector"),
          function(x, v, t = 1.0, u = NULL,
                   Markov = FALSE, transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@x,
                              ia = x@i + 1,
                              ja = x@p + 1,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'CCS',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

##' @rdname expv-methods
##' @aliases expv,TsparseMatrix,vector-method
##' @importClassesFrom Matrix TsparseMatrix
setMethod("expv", signature("TsparseMatrix", "vector"),
          function(x, v, t = 1.0, u = NULL,
                   Markov = FALSE, transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@x,
                              ia = x@i + 1,
                              ja = x@j + 1,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'COO',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

##' @rdname expv-methods
##' @aliases expv,matrix.csc,vector-method
##' @importClassesFrom SparseM matrix.csc
##' @importMethodsFrom SparseM t
setMethod("expv", signature("matrix.csc", "vector"),
          function(x, v, t = 1.0, u = NULL,
                   Markov = FALSE, transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@ra,
                              ia = x@ja,
                              ja = x@ia,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'CCS',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

##' @rdname expv-methods
##' @aliases expv,matrix.csr,vector-method
##' @importClassesFrom SparseM matrix.csr
setMethod("expv", signature("matrix.csr", "vector"),
          function(x, v, t = 1.0, u = NULL,
                   Markov = FALSE, transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@ra,
                              ia = x@ia,
                              ja = x@ja,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'CRS',
                              ...)
              v <- w[, k]
            }
            w
          }
          )

##' @rdname expv-methods
##' @aliases expv,matrix.coo,vector-method
##' @importClassesFrom SparseM matrix.coo
setMethod("expv", signature("matrix.coo", "vector"),
          function(x, v, t = 1.0, u = NULL,
                   Markov = FALSE, transpose = Markov, ...) {
            if (transpose)
              x <- t(x)
            w <- matrix(0, length(v), length(t))
            if (length(t) > 1)
              t <- c(t[1], diff(t))
            for(k in seq_along(t)) {
              w[, k] <- Rexpv(a = x@ra,
                              ia = x@ia,
                              ja = x@ja,
                              n = nrow(x),
                              v = v,
                              u = u,
                              t = t[k],
                              Markov = Markov,
                              storage = 'COO',
                              ...)
              v <- w[, k]
            }
            w
          }
          )




            
