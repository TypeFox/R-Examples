# @title Log-likelihood function of a CUBE model without covariates
# @aliases loglikcube
# @description Compute the log-likelihood function of a CUBE model without covariates fitting 
# the given absolute frequency distribution.
# @usage loglikcube(m,freq,pai,csi,phi)
# @param m Number of ordinal categories
# @param freq Vector of the absolute frequency distribution
# @param pai Uncertainty parameter
# @param csi Feeling parameter
# @param phi Overdispersion parameter
#' @keywords internal


loglikcube <-
function(m,freq,pai,csi,phi){t(freq)%*%log(probcube(m,pai,csi,phi))}
