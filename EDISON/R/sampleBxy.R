#' Sample regression coefficients.
#' 
#' This function samples the regression coefficients given the current state of
#' the MCMC simulation.
#' 
#' 
#' @param xi Response data.
#' @param y Target data.
#' @param Sig2 Sigma squared.
#' @param delta2 Signal-to-noise hyperparameter.
#' @return The regression parameters.
#' @author Sophie Lebre
#' @references For details of the regression model, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export sampleBxy
sampleBxy <-
function(xi, y, Sig2, delta2){
  # INPUT: xi, yi, sig2i, delta2
  # OUTPUT: B
  Ml = (delta2 / (delta2+1)) * pseudoinverse(t(xi) %*% xi)
  out = mvrnorm(1, mu=Ml %*% t(xi) %*% y, Sigma=Sig2*Ml)
  return(out)
}

