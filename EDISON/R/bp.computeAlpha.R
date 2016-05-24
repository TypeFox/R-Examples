#' Computes the acceptance ratio of two changepoint configurations.
#' 
#' This function computes the acceptance ratio of two changepoint
#' configurations with networks in a changepoint birth or death move.
#' 
#' 
#' @param birth \code{1} for a changepoint birth move, \code{-1} for a
#' changepoint death move.
#' @param lNew Number of edges in the new segment.
#' @param kminus Minimal number of changepoints between the two compared models
#' (equal to \code{s} for a birth move, \code{s-1} for a death move.
#' @param Ekl Changepoint on the left of proposed changepoint.
#' @param Estar Changepoint being inserted or deleted.
#' @param Ekr Changepoint on the right of proposed changepoint.
#' @param yL Response data (left).
#' @param PxL Projection matrix (left).
#' @param yR Response data (right).
#' @param PxR Projection matrix (right).
#' @param y2 Response data (both).
#' @param Px2 Projection matrix (both).
#' @param D Hyperparameters for the number of edges in each segment.
#' @param delta2 Hyperparameters for the empirical covariance (signal-to-noise
#' ratio).
#' @param q Total number of nodes in the network.
#' @param smax Maximum number of changepoints.
#' @param v0 Hyperparameter.
#' @param gamma0 Hyperparameter.
#' @param prior_ratio Ratio of network structure priors.
#' @author Sophie Lebre
#' @seealso \code{\link{cp.birth}}, \code{\link{cp.death}}
#' @references For more information about the model, see:
#' 
#' Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian networks with
#' Bayesian regularization for inferring gene regulatory networks with
#' gradually time-varying structure", Machine Learning.
#' @export bp.computeAlpha
bp.computeAlpha <-
function(birth,lNew,kminus,Ekl,Estar,Ekr,yL,PxL,yR,PxR,y2,Px2,D,delta2, q, smax, v0, gamma0, prior_ratio=1){
  # birth = 1 for birth, -1 for death.
  # lNew = number of edges in the new phase
  # kminus = minimal number of changepoints between the 2 compared models (=s for birth, s-1 for death) -> INUTILE !!!
  # Ekl = 
  # Estar =
  # yL,yR, y2 : response data  (left, right, both)
  # PxL, PxR, Px2 : projection matrix (left, right, both)
  # D : hyperparms for the number of edges in each phase. (Number of edges s ~ truncated Poisson P(D). )
  # delta2 : hyperparms for empirical covariance (can be seen as the expected  signal-to-noise ratio)  ~IG(alphad2,betad2)

  logR=  +(v0/2)*log(gamma0/2)-lgamma(v0/2)  -log( (sqrt(delta2+1))^(lNew+1) ) +lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-lgamma(((Ekr-Ekl)+v0)/2) +(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))  +(((Ekr-Estar))/2)*log (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))  -(v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)
  
  logR = logR + log(prior_ratio)
  
  logR = birth*logR
  
  if(logR > 0){
    res=1
  } else {
    res=exp(logR)
  }
  
  return(res)
}

