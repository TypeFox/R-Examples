# @title Log-likelihood function for an IHG model without covariates
# @aliases loglikihg
# @description Compute the log-likelihood function for an IHG model without covariates for
#  the given absolute frequency distribution.
#  @usage loglikihg(m, freq, theta)
# @param m Number of ordinal categories
# @param freq Vector of the absolute frequency distribution
# @param theta Preference parameter
# @seealso loglikIHG
#' @keywords internal




loglikihg <-
function(m,freq,theta){
  t(freq)%*%log(probihg(m,theta))
}
