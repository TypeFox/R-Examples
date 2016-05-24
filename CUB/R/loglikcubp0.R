# @title Log-likelihood function of a CUB model with covariates for the uncertainty component
# @description Compute the log-likelihood function of a CUB model fitting ordinal responses with covariates 
# for explaining the uncertainty component.
# @aliases loglikcubp0
# @usage loglikcubp0(m, ordinal, Y, bbet, ccsi)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of selected covariates for explaining the uncertainty component
# @param bbet Vector of parameters for the uncertainty component, with length equal to 
# NCOL(Y)+1 to account for an intercept term (first entry of bbet)
# @param ccsi Feeling parameter
#' @keywords internal


loglikcubp0 <-
function(m,ordinal,Y,bbet,ccsi){
  prob<-probbit(m,ccsi)
  probn<-prob[ordinal]
  eta<-logis(Y,bbet)
  return(sum(log(eta*(probn-1/m)+1/m)))
}
