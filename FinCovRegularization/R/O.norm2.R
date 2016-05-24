#' @title The Squared Operator Norm
#'
#' @description
#' Calculate the squared Operator norm of a matrix
#'
#' @param matrix a matrix
#' @return a scalar of the squared Operator norm
#' @examples
#' data(m.excess.c10sp9003)
#' cov.SAM <- cov(m.excess.c10sp9003)
#' O.norm2(cov.SAM)
#' @export

O.norm2 <- function(matrix) {
  O2 <- norm(matrix,"2")^2
  return(O2)
}