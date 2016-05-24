#' Orthogonal Procrustes Analysis
#' 
#' Simple function to rotate matrix X so that it matches the target matrix Z as
#' closely as possible, by minimizing ||Z - XQ|| where Z and X are of the same
#' size and Q is an orthogonal matrix. The algorithm is based on the singular
#' value decomposition (SVD) (see e.g. the reference).
#' 
#' @param Z The target matrix
#' @param X The matrix to be rotated, which must be of the same size as Z.
#' @return A list with the following 2 elements: \item{Q}{The rotation matrix}
#' \item{XQ}{The matrix X after rotation} 
#' @references Gower, J. C. and Hand, D.J. (1996). Biplots (Vol. 54). CRC Press.
#' @export orthprocr
orthprocr <- function(Z, X){
  stopifnot(all.equal(dim(Z), dim(X)))
  ZX <- crossprod(Z, X)
  svdZX <- svd(ZX)
  Q <- tcrossprod(svdZX$v, svdZX$u)
  list(Q = Q, XQ = X %*% Q)
}
