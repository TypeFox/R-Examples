#' Computing Wasserstein Metric
#' @param cases name of case group data (matrix sample * feature).
#' @param control names of control group data (matrix sample * feature).
#' @param paranum the number of quatile discretization + 1. Default is discretized by 1\%. 
#' @param q power of Wasserstein metric. Default is q = 2.
#' @return vector of Wasserstein metric
#' @author Yusuke Matsui & Teppei Shimamura
#' @examples
#' nrep <- 12
#' cases <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep))
#' cases <- do.call("rbind",cases)
#' control <- Map(rbeta,rep(30,nrep),rep(1,nrep),rep(5,nrep))
#' control <- do.call("rbind",control)
#' wasserMetric(cases,control)
#' 
#' 
#' @export
#' 

wasserMetric <- function(cases,control,paranum = 101, q = 2){
  
  d <- wasserCpp_mat(cases,control,paranum = paranum, q = q)
  
  d
}
