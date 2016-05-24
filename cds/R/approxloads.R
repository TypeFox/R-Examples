#' Low Rank Approximation LL' of a Square Symmetrix Matrix R
#' 
#' Uses the eigendecomposition of a square, symmetrix matrix R to obtain the
#' loadings matrix L such that R is approximated by LL', with L restricted to
#' have \code{r} columns. Hence LL' is a rank \code{r} approximation of R. The
#' eigendecomposition of R is used to obtain L from the first \code{r}
#' eigenvectors and eigenvalues. In case \code{procr.target} is not
#' \code{NULL}, L is further rotated through orthogonal Procrustes analysis to
#' match as closely as possible the matrix \code{procr.target} through
#' \code{\link{orthprocr}}.
#' 
#' @param R Square, symmetric matrix R to be approximated
#' @param r The required rank of the approximation
#' @param procr.target Optional; the target matrix for L in the orthogonal
#' Procrustes analysis
#' @param refl.target Optional; the matrix to check against for possible
#' reflections of the loading vectors.
#' @examples
#' R <- rcormat(10, r = 3)
#' all.equal(R$L, approxloads(R$R, r = 3, procr.target = R$L))
#' @export approxloads
approxloads <- function(R, r = 3, procr.target = NULL, refl.target = NULL){
  stopifnot(is.matrix(R), ncol(R) == nrow(R), all.equal(R, t(R)))
  if(any(is.na(R))) return(NA)
  eigR <- eigen(R)
  L <- eigR$vectors[,1:r, drop = FALSE] %*% (sqrt(eigR$values[1:r])*diag(, nrow = r, ncol = r))
  if(!is.null(procr.target)){
    L <- orthprocr(Z = procr.target, X = L)$XQ
  } 
  if(!is.null(refl.target)){
      r.refl <- min(r, ncol(refl.target))
      RMSEorig <- colSums((L[, 1:r.refl, drop = FALSE] - refl.target[, 1:r.refl, drop = FALSE])^2)
      RMSErefl <- colSums((L[, 1:r.refl, drop = FALSE] + refl.target[, 1:r.refl, drop = FALSE])^2)
      ind <- apply(cbind(RMSEorig, RMSErefl), 1, which.min)
      L[, 1:r.refl][, ind == 2] <- -1*L[, 1:r.refl][, ind == 2]
  }
  L
}
