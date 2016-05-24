# @title Log-likelihood function of a CUB model without covariates
# @description Compute the log-likelihood function of a CUB model without covariates for a given 
# absolute frequency distribution.
# @aliases loglikcub00
# @usage loglikcub00(m, freq, pai, csi)
# @param m Number of ordinal categories
# @param freq Vector of the absolute frequency distribution
# @param pai Uncertainty parameter
# @param csi Feeling parameter
#' @keywords internal


loglikcub00 <-
function(m,freq,pai,csi){t(freq)%*%log(probcub00(m,pai,csi))}
