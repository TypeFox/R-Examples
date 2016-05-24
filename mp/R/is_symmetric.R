#' Tests whether the given matrix is symmetric.
#'
#' @param mat Matrix to be tested for symmetry.
#' @return Whether the matrix is symmetric.
#'
#' @export
is.symmetric <- function(mat) {
  if (!is.matrix(mat)) {
    stop("supplied argument is not a matrix")
  }

  if (nrow(mat) != ncol(mat)) {
    stop("given matrix is not square")
  }

  tmat <- t(mat)
  identical(mat[lower.tri(mat)], tmat[lower.tri(mat)])
}
