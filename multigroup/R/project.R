#' @title Project operatior
#' @description 
#' To calculate projection of a matrix
#' @param x a numeric vector
#' @export
#' @keywords internal
project <- function(X) {
  # projection
  Y = X %*% ginv(t(X) %*% X) %*% t(X)
}

