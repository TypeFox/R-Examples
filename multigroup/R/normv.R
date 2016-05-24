#' @title Normalize function for a vector
#' @description 
#' To calculate normalize of a vector
#' @param x a numeric vector
#' @return a vector with norm 1
#' @export
#' @keywords internal
normv <- function(x){
  
  normx = sqrt(sum(x * x))
  y    = x / normx 
  return(y)
}
