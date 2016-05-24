# @title Log-likelihood function of a CUBE model with covariates
# @aliases loglikcubecov
# @description Compute the log-likelihood function of a CUBE model for ordinal responses,
#  with covariates for explaining all the three parameters.
#  @usage loglikcubecov(m, ordinal, Y, W, Z, bet, gama, alpha)
# @param m Number of ordinal categories
# @param ordinal  Vector of ordinal responses
# @param Y Matrix of covariates for explaining the uncertainty component
# @param W Matrix of covariates for explaining the feeling component
# @param Z Matrix of covariates for explaining the overdispersion component
# @param bet Vector of parameters for the uncertainty component, with length equal to 
# NCOL(Y) + 1 to account for an intercept term (first entry of bet)
# @param gama Vector of parameters for the feeling component, with length equal to 
# NCOL(W) + 1 to account for an intercept term (first entry of gama)  
# @param alpha Vector of parameters for the overdispersion component, with length equal to 
# NCOL(Z) + 1 to account for an intercept term (first entry of alpha) 
#' @keywords internal



loglikcubecov <-
function(m,ordinal,Y,W,Z,bet,gama,alpha){
  paivett<-logis(Y,bet); csivett<-logis(W,gama); 
  phivett<-1/(-1+ 1/(logis(Z,alpha))) 
  probi<-paivett*(betabinomial(m,ordinal,csivett,phivett)-1/m)+1/m
  return(sum(log(probi)))
}
