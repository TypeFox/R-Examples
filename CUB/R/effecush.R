# @title Auxiliary function for the log-likelihood estimation of CUSH models with covariates
# @description Compute the opposite of the loglikelihood function for CUSH models
# with covariates to explain the shelter effect.
# @aliases effecush
# @usage effecush(paravec, esternocush, shelter, m)
# @param paravec Vector of the initial parameters estimates
# @param esternocush Matrix binding together the vector of ordinal data and the matrix XX of explanatory
# variables whose first column is a column of ones needed to consider an intercept term
# @param shelter Category corresponding to the shelter choice
# @param m Number of ordinal categories
# @details It is called as an argument for "optim" within CUSH function (when no covariate is included)
#  as the function to minimize.
#' @keywords internal 


effecush <-
function(paravec,esternocush,shelter,m){
  ordinal<-esternocush[,1]
  ncovar<-ncol(esternocush)
  X<-esternocush[,3:ncovar]
  return(-loglikcushcov(m,ordinal,X,paravec,shelter))
}
