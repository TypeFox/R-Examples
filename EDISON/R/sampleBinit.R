#' Sample initial regression coefficients.
#' 
#' This function samples the initial regression coefficients for the networks.
#' 
#' 
#' @param Si Network structure.
#' @param sig2 Sigma squared.
#' @param delta2 Signal-to-noise ratio hyperparameter.
#' @param X Input data.
#' @param q Number of nodes.
#' @return Returns a vector of regression coefficients.
#' @author Sophie Lebre
#' @references For details of the regression model, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export sampleBinit
sampleBinit <-
function(Si, sig2, delta2, X, q){
  ### INPUT: Si=S[i,],sig2=Sig2[i],delta2,
  ###        X the observed data for predictors.
  ###        q number of predictors
  ### OUTPUT: vector newB.
  ### depends on: q the number of predictors.
  newB <- array(0,q+1)
  
  for(l in which(Si == 1)){
    newB[l] <- rnorm(1, mean=0, sd=sqrt(delta2 * sig2 * t(X[,l]) %*% X[,l]))
  }
  
  return(newB)
}

