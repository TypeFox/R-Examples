#' @title Computes quadratic form x' A x
#' 
#' @description
#' \code{quadratic_form} computes the quadratic form \eqn{\mathbf{x}' \mathbf{A} \mathbf{x}} for an
#' \eqn{n \times n} matrix \eqn{\mathbf{A}} and an \eqn{n}-dimensional vector
#' \eqn{\mathbf{x}}, i.e., a wrapper for \code{t(x) \%*\% A \%*\% x}. 
#' 
#' \code{fill_symmetric} and \code{quadratic_form} work with 
#' real and complex valued matrices/vectors.
#' 
#' @param mat numeric; \eqn{n \times n} matrix (real or complex).
#' @param vec numeric; \eqn{n \times 1} vector (real or complex).
#' @return 
#' A real/complex value \eqn{\mathbf{x}' \mathbf{A} \mathbf{x}}.
#' @export
#' @keywords math univar
#' @examples
#'  set.seed(1)
#'  AA <- matrix(1:4, ncol = 2)
#'  bb <- matrix(rnorm(2))
#'  t(bb) %*% AA %*% bb
#'  quadratic_form(AA, bb)
#' 
#' 

quadratic_form <- function(mat, vec) {
  # computes the quadratic form vec' * mat * vec
  
  stopifnot(is.matrix(mat))
  dim.mat <- dim(mat)
  stopifnot(dim.mat[1] == dim.mat[2])
  
  vec <- matrix(vec)
  qp <- crossprod(Conj(vec), crossprod(mat, vec))[1, 1]
  # convert to real if imaginary part is 0
  if (round(Im(qp), 6) == 0) {
    qp <- Re(qp)
  }
  names(qp) <- NULL
  return(qp)
} 

#' @rdname quadratic_form
#' @export
#' @description
#' \code{fill_hermitian} fills up the lower triangular part (\code{NA})
#' of an upper triangular matrix to its
#' Hermitian (symmetric if real matrix) version, such that it satisfies 
#' \eqn{\mathbf{A} = \bar{\mathbf{A}}'}, where \eqn{\bar{z}} is the complex
#' conjugate of \eqn{z}.  If the matrix is real-valued this makes it 
#' simply symmetric.
#' 
#' Note that the input matrix must have a \strong{real-valued} diagonal and 
#' \code{NA}s in the lower triangular part.
#' @examples
#' 
#' AA <- matrix(1:16, ncol = 4)
#' AA[lower.tri(AA)] <- NA
#' AA
#' 
#' fill_hermitian(AA)
#' 

fill_hermitian <- function(mat) {
  stopifnot(inherits(mat, "matrix"),
            ncol(mat) == nrow(mat),
            all.equal(Im(diag(mat)), rep(0, ncol(mat))))
  
  if (!identical(dim(mat), c(1, 1))) {
    lower.ind <- lower.tri(mat)
    mat.lower <- mat[lower.ind]
    # lower triangular part must be NA
    stopifnot(all(is.na(mat.lower)))
  
    mat[lower.ind] <- Conj(t(mat)[lower.ind])
  }
  return(mat)
} 