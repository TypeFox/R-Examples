#' Update sigma squared variances.
#' 
#' This function samples new values for the sigma squared variances, given the
#' current network structure. A multivariate distribution is assumed.
#' 
#' 
#' @param phase Current segment.
#' @param X Input response data.
#' @param Y Input target data.
#' @param E Changepoints.
#' @param Sall Network structure.
#' @param Ball Regression coefficients.
#' @param Sig2 Current sigma squared values.
#' @param Mphase Segment positions.
#' @param alphad2 Hyperparameter for gamma prior.
#' @param betad2 Hyperparameter for gamma prior.
#' @param v0 Hyperparameter for inverse gamma prior.
#' @param gamma0 Hyperparameter for inverse gamma prior.
#' @return The new samples sigma squared values.
#' @author Sophie Lebre
#' @seealso \code{\link{updateSigSolo}}
#' @references For more information about the model, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export updateSigMulti
updateSigMulti <-
function(phase, X, Y, E, Sall, Ball, Sig2, 
                           Mphase, alphad2, betad2, v0, gamma0){
  
  posPhase = which(E == phase)
  S = Sall[posPhase,]
  k = sum(Sall[posPhase,])-1
  
  y = Y[(Mphase[phase]:(Mphase[E[posPhase+1]]-1))]
  x = X[(Mphase[phase]:(Mphase[E[posPhase+1]]-1)),]
  delta2 = rinvgamma(1, shape=k + alphad2, scale=betad2 + Ball[posPhase, which(S == 1)] %*% t(x[, which(S == 1)]) %*% x[,which(S == 1)] %*% Ball[posPhase,which(S == 1)] / (2 * Sig2) )
  matPx = computePx(length(y), x[, which(S == 1)], delta2)
    
  total = t(y) %*% matPx %*%y
  out = rinvgamma(1, shape=v0/2 + length(y)/2, scale=(gamma0 + total)/2)
  return(out )
  
}

