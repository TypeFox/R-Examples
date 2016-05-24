#' Quadratic form of a matrix and a vector
#'
#' We compute the quadratic form of a vector and a matrix in an efficient
#' manner. Let \code{x} be a real vector of length \code{p}, and let \code{A} be
#' a p x p real matrix. Then, we compute the quadratic form \eqn{q = x' A x}.
#'
#' A naive way to compute the quadratic form is to explicitly write
#' \code{t(x) \%*\% A \%*\% x}, but for large \code{p}, this operation is
#' inefficient. We provide a more efficient method below.
#'
#' Note that we have adapted the code from:
#' \url{http://tolstoy.newcastle.edu.au/R/help/05/11/14989.html}
#' 
#' @param A matrix of dimension p x p
#' @param x vector of length p
#' @return scalar value
quadform <- function(A, x) {
  drop(crossprod(x, A %*% x))
}

#' Quadratic Form of the inverse of a matrix and a vector
#'
#' We compute the quadratic form of a vector and the inverse of a matrix in an
#' efficient manner. Let \code{x} be a real vector of length \code{p}, and let
#' \code{A} be a p x p nonsingular matrix. Then, we compute the quadratic form
#' \eqn{q = x' A^{-1} x}.
#'
#' A naive way to compute the quadratic form is to explicitly write
#' \code{t(x) \%*\% solve(A) \%*\% x}, but for large \code{p}, this operation is
#' inefficient. We provide a more efficient method below.
#'
#' Note that we have adapted the code from:
#' \url{http://tolstoy.newcastle.edu.au/R/help/05/11/14989.html}
#' 
#' @param A matrix that is p x p and nonsingular
#' @param x vector of length p
#' @return scalar value
quadform_inv <- function(A, x) {
  drop(crossprod(x, solve(A, x)))
}

#' Centers the observations in a matrix by their respective class sample means
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @return matrix with observations centered by its corresponding class sample
#' mean
center_data <- function(x, y) {
  x <- as.matrix(x)
  y <- as.factor(y)

  # Notice that the resulting centered data are sorted by class and do not
  # preserve the original ordering of the data.
  x_centered <- tapply(seq_along(y), y, function(i) {
    scale(x[i, ], center = TRUE, scale = FALSE)
  })
  x_centered <- do.call(rbind, x_centered)

  # Sorts the centered data to preserve the original ordering.
  orig_ordering <- do.call(c, tapply(seq_along(y), y, identity))
  x_centered[orig_ordering, ] <- x_centered
  x_centered
}

#' Computes the inverse of a symmetric, positive-definite matrix using the
#' Cholesky decomposition
#'
#' This often faster than \code{\link{solve}} for larger matrices.
#' See, for example:
#' \url{http://blog.phytools.org/2012/12/faster-inversion-of-square-symmetric.html}
#' and
#' \url{http://stats.stackexchange.com/questions/14951/efficient-calculation-of-matrix-inverse-in-r}.
#'
#' @export 
#' @param x symmetric, positive-definite matrix
#' @return the inverse of \code{x}
solve_chol <- function(x) {
  chol2inv(chol(x))
}
