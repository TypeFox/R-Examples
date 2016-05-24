#'  Generate a random sample from a truncated normal distribution.
#'  
#' @param mean the mean of the normal distribution
#' @param variance The variance of the normal distribution
#'    
#' @return A random sample from the specified truncated normal distribution
#'  
#' @example tests/testthat/examples_fcn_doc/examples_getTruncatedNormal.R
#' @export
#' @keywords internal
getTruncatedNormal <- function(mean,variance){
  while(TRUE){
    n = mean+randn(1,1)*sqrt(variance)
    if((sign(n)==sign(mean))){
      break
    }
  }
  return( n) 
}
