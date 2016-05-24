#'Power for matrix elements
#'
#'Set each element of matrix \code{x} to the power \code{n}
#'
#'@param x a matrix
#'
#'@param n a real number
#'
#'@keywords internal

"%^%" <- function(x, n){
  with(eigen(x), vectors %*% (values^n * t(vectors)))
}
