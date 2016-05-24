#' Wasserstein Metric
#' @param cases name of case group data (matrix sample * feature).
#' @param control names of control group data (matrix sample * feature).
#' @param paranum the number of quatile discretization + 1. Default is discretized by 1\%. 
#' @param q power of Wasserstein metric. Default is q = 2.
#' @return Wasserstein metric
#' @author Yusuke Matsui & Teppei Shimamura
#' @examples
#' cases <- rbeta(30,1,5)
#' control <- rbeta(30,2,5)
#' wasserMetric.v(cases,control)
#' 
#' @export
#' 

wasserMetric.v <- function(cases,control,paranum = 101, q = 2){
  #library(Rcpp)
  
  d <- wasserCpp(cases,control,paranum = paranum, q = q)
  
  d
}
