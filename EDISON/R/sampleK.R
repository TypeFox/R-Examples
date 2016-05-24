#' Sample initial number of changepoints.
#' 
#' This function samples the initial number of changepoints from a sparse
#' Poisson prior.
#' 
#' 
#' @param mini Minimum value.
#' @param maxi Maximum value.
#' @param lambda Parameter of the Poisson distribution.
#' @param nb Number of values to sample.
#' @return The sampled number of changepoints.
#' @author Sophie Lebre
#' @references For more information on the prior choice and sampling, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export sampleK
sampleK <-
function(mini, maxi, lambda, nb){
  if( mini == maxi) { 
    print("Error with sampling from a truncated Poisson: min value = max value.") }
  
  out = sample(mini:maxi, nb, replace=TRUE, 
               prob=lambda^(mini:maxi)/apply(matrix(mini:maxi, 1, maxi-mini+1),
                                             2, factorial))
  return(out)
}

