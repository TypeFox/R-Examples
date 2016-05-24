# =========== Reliability estimators ===========

#'@title Reliabilities as products of weights and loadings
#'
#'@description
#'Calculates reliabilities as a matrix product of loadings and weights.
#'
#'@param S the data covariance matrix
#'
#'@param loadings matrix of factor loading estimates
#'
#'@param W matrix of weights
#'
#'@param ... All other arguments are ignored.
#'
#'@return a named vector of estimated composite reliabilities.
#'
#'@family reliability estimators
#'
#'@export

reliability.weightLoadingProduct <- function(S, loadings, W, ...){
  
  Q <- diag(W %*% loadings)^2
  
  # Any composite with no reflective indicators is fixed to be perfectly reliable
  Q[apply(loadings==0,2,all)] <- 1
  
  Q
}
