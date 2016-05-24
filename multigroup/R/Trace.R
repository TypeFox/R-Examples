#' @title trace of a matrix
#' @description 
#' To calculate trace
#' @param x a numeric vector
#' @export
#' @keywords internal
Trace <- function(X){
  
  Y = sum( diag( X %*% t(X)))
}
