#' @title Normalize function for a matrix
#' 
#' @description 
#' To calculate normmalize of a matrix
#' 
#' @param X a numeric matrix
#' @return a matrix with norm 1
#' @export
#' @keywords internal
normM <- function(X){
  
  if (class(X) == 'data.frame') {
    X = as.matrix(X)
  }
  
  normX = sqrt(sum(diag(t(X) %*% X)))
  Y     = X / normX 
  return(Y)
}