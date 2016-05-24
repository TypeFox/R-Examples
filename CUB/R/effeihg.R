# @title Auxiliary function for the log-likelihood estimation of IHG models without covariates
# @description Compute the opposite of the log-likelihood function for an IHG model without covariates.  
# @aliases effeihg
# @usage effeihg(theta, m, freq)
# @param theta Initial estimate for the parameter of the IHG distribution
# @param m Number of ordinal categories
# @param freq Vector of the absolute frequency distribution of the ordinal responses
# @details It is called as an argument for "optim" within IHG function (when no covariate is specified)
# as the function to minimize.
#' @keywords internal 


effeihg <-
function(theta,m,freq){-loglikihg(m,freq,theta)}
