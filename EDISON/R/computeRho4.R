#' Calculate proposal frequencies for changepoint moves.
#' 
#' This function calculates the frequency at which each of the different
#' changepoint moves is proposed. For the poisson network structure prior, this
#' ensures that the proposal frequency is equal to the prior probability.
#' 
#' 
#' @param k The number of hidden states.
#' @param kmin Minimum number of hidden states.
#' @param kmax Maximum number of hidden states
#' @param c Parameter.
#' @param lambda Hyperparameter controlling the number of hidden states.
#' @return Vector containing the proposal frequencies for the different
#' changepoint moves.
#' @author Sophie Lebre
#' @references For more information about the hyperparameters and the
#' functional form of the likelihood, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export computeRho4
computeRho4 <-
function(k, kmin, kmax, c, lambda){
  # INPUT:   k the number of hidden states
  #          kmax the maximal number of hidden states
  #          c constante
  #          lambda simulation parameter for the number of hidden states
  # OUTPUT: rho

  if(c == 0) return(c(0, 0, 0))
  
  rho = array(1,4)
  
  if(k == kmax) { rho[1] = 0 } else { rho[1] = c*min(1,lambda/(k+1)) }
  if(k == kmin) { rho[2] = rho[1] } else { rho[2] = rho[1]+c*min(1,k/lambda) }
  if(k > 0) { rho[3] = rho[2]+(1-rho[2])/3 } else{ rho[3] = rho[2] }

  return(rho)
}

