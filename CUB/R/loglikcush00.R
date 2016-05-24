# @title Log-likelihood function for a CUSH model without covariates
# @aliases loglikcush00
# @description Compute the log-likelihood function for a CUSH model 
# without covariate for the given ordinal responses.
# @usage loglikcush00(m,ordinal,delta,shelter)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param delta Shelter parameter
# @param shelter Category corresponding to the shelter choice
# @seealso \code{\link{CUSH}}
#' @keywords internal


loglikcush00<-function(m,ordinal,delta,shelter){
  n<-length(ordinal);
  freq<-tabulate(ordinal,nbins=m)
  fc<-freq[shelter]/n
  loglik<-n*((1-fc)*log(1-delta)+fc*log(1+(m-1)*delta)-log(m))
  return(loglik)
}

