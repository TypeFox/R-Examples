# @title Log-likelihood function of a CUB model with covariates for both feeling and uncertainty
# @description Compute the log-likelihood function of a CUB model fitting ordinal data
#  with covariates for explaining both the feeling and the uncertainty components.
# @aliases loglikcubpq
# @usage loglikcubpq(m, ordinal, Y, W, bet, gama)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of selected covariates for explaining the uncertainty component
# @param W Matrix of selected covariates for explaining the feeling component
# @param bet Vector of parameters for the uncertainty component, with length equal to 
# NCOL(Y)+1 to account for an intercept term (first entry of bbet)
# @param gama Vector of parameters for the feeling component, whose length equals
# NCOL(W) + 1 to account for an intercept term (first entry of gama)
#' @keywords internal


loglikcubpq <-
function(m,ordinal,Y,W,bet,gama){
  probn<-probcubpq(m,ordinal,Y,W,bet,gama)
  return(sum(log(probn)))
}
