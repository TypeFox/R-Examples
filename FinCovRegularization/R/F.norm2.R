#' @title The Squared Frobenius Norm
#'
#' @description
#' Calculate the squared Frobenius norm of a matrix
#'
#' @param matrix a matrix
#' @return a scalar of the squared Frobenius norm
#' @examples
#' data(m.excess.c10sp9003)
#' cov.SAM <- cov(m.excess.c10sp9003)
#' F.norm2(cov.SAM)
#' @export

F.norm2 <- function(matrix) {
  F2 <- norm(matrix,"F")^2
  return(F2)
}