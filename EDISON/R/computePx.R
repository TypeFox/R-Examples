#' Compute projection matrix.
#' 
#' This function computes the projection matrix that is needed for calculation
#' of the likelihood.
#' 
#' 
#' @param len Delimiting breakpoints.
#' @param x The observations of x in the corresponding state.
#' @param delta2 Signal-to-noise ratio hyperparameter.
#' @return The projection matrix.
#' @author Sophie Lebre
#' @seealso \code{\link{CalculateLikelihoodRatio}}
#' @references For more information about the hyperparameters and the
#' functional form of the likelihood, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export computePx
computePx <-
function(len, x, delta2){
  # INPUT: len, delimiting breakpoints.
  #        x, the observations of X in the corresponding state
  #        delta2.
  # OUTPUT: the projection matrix Px.

  if(prod(dim(x))>0) {
    moins = (delta2/(delta2+1))* x%*%pseudoinverse(t(x)%*%x)%*%t(x)
  } else {
    moins = matrix(0,len,len)
  }
  
  Px = diag(1,len)-moins
  
  return(Px)
}

