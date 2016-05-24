#' Logistic detection function
#'
# Logistic detection function with x = distance and z = other covariates,
# X(x,z) is design matrix for x,z.
#
# h(x,z) = exp (X(x,z)*beta)/(1+exp(X(x,z)*beta)
# g(x,z)= h(x,z)/h(0,z)
#'
#' @param distance perpendicular distance vector
#' @param theta scale parameters
#' @param w scale covariate matrix
#' @param std if TRUE uses scale=1
#'
#' The routine returns a vector of probabilities that the observation
#' were detected given they were at the specified distance and assuming that
#' g(0)=1 (ie a standard line transect detection function).
logisticdetfct <- function(distance, theta, w, std=FALSE){
  exp(as.matrix(w[[1]]) %*% beta)/
    (1 + exp(as.matrix(w[[1]]) %*% beta))/
    g0(beta=theta,w[[2]])
}
