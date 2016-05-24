#' @title Deflation
#' 
#' @description 
#' To deflate a matrix
#' 
#' @param X a numeric matrix
#' @param d a vector to deflate
#' @return Y a deflated matrix
#' @export
#' @keywords internal
deflation <- function(X,d){
  Y = X -   d %*% t(d) %*% X / sum(d * d)
  return(Y) 
   
}