# @title Log-likelihood function for a CUSH model with covariates
# @aliases loglikcushcov
# @description Compute the log-likelihood function for a CUSH model 
# with covariates for the given ordinal responses.
# @usage loglikcushcov(m, ordinal, X, omega, shelter)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param X Matrix of selected covariates for explaining the shelter parameter
# @param omega Vector of parameters for explaining the shelter effect, 
# with length equal to NCOL(X)+1 to account for an intercept term (first entry of omega)
# @param shelter Category corresponding to the shelter choice
# @seealso \code{\link{CUSH}}
#' @keywords internal



loglikcushcov <-
function(m,ordinal,X,omega,shelter){
  deltavett<-logis(X,omega)
  dummy<-ifelse(ordinal==shelter,1,0)
  probi<-deltavett*(dummy-1/m)+1/m
  return(sum(log(probi)))
}
