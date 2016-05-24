#' Sample delta squared.
#' 
#' This function samples the signal-to-noise hyperparameter delta squared.
#' 
#' 
#' @param pos The current segment.
#' @param x Data,
#' @param q Number of nodes.
#' @param B Regression coefficients.
#' @param S Network structure.
#' @param sig2 Sigma squared.
#' @param alphad2 Gamma prior hyperparameter.
#' @param betad2 Gamma prior hyperparameter.
#' @return New sample of delta squared.
#' @author Sophie Lebre
#' @references For details of the sampling, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export sampleDelta2
sampleDelta2 <-
function(pos, x, q, B, S, sig2, alphad2, betad2){
  # INPUT: pos, the considered state
  #        xPos, the observations of X in state i
  #        B,S,Sig2
  #        alphad2,betad2.
  # OUTPUT: delta2
 
  plus = 0
  
  if(sum(S[pos,]) > 0){
    Bi = B[pos, which(S[pos,] == 1)]
    xi = x[, which(S[pos,] == 1)]
    plus = Bi %*% t(xi) %*% xi %*% Bi / (2* sig2)
  }
  
  out = rinvgamma(1, shape=sum(S[pos,1:q]) + alphad2, 
                  scale=betad2 + plus)
  
  return(out)
}

