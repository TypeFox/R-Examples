# @title Auxiliary function for the log-likelihood estimation of IHG models with covariates
# @description Compute the opposite of the log-likelihood function for an IHG model with covariates
# for the preference parameter.  
# @aliases effeihgcov
# @usage effeihgcov(nu, ordinal, U, m)
# @param  nu Vector of the starting values for the parameters to be estimated, with length equal to  
#  NCOL(U)+1 to account for an intercept term (first entry of \eqn{nu})
# @param ordinal Vector of ordinal responses
# @param U Matrix of the explanatory variables for the preference parameter \eqn{\theta}
# @param m Number of ordinal categories
# @details It is called as an argument for "optim" within IHG function (with covariates)
#  as the function to minimize.
#' @keywords internal 

effeihgcov <-
function(nu,ordinal,U,m){
  -loglikihgcov(m,ordinal,U,nu)
}
