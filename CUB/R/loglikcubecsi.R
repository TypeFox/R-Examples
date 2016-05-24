# @title Log-likelihood function of CUBE model with covariates only for feeling
# @aliases loglikcubecsi
# @description 
# Compute the log-likelihood function of a CUBE model for ordinal data with subjects' 
# covariates only for feeling.
# @usage loglikcubecsi(m, ordinal, W, pai, gama, phi)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param W Matrix of covariates for explaining the feeling component
# @param pai Uncertainty parameter
# @param gama Vector of parameters for the feeling component, with length 
# equal to NCOL(W) + 1 to account for an intercept term (first entry of gama)
# @param phi Overdispersion parameter
# @seealso loglikCUBE
#' @keywords internal


loglikcubecsi <-
function(m,ordinal,W,pai,gama,phi){
  csivett<-logis(W,gama)
  probi<-pai*(betabinomialcsi(m,ordinal,csivett,phi)-1/m)+1/m
  return(sum(log(probi)))
}
