#' Compute value of p(0) using a logit formulation
#'
#' @param beta logistic parameters
#' @param z design matrix of covariate values
#'
#' @return vector of p(0) values
#' @author Jeff Laake
g0 <- function(beta, z){
  exp(z %*% beta)/(1 + exp(z %*% beta))
}
