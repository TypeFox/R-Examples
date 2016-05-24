#' -2 * log-likelihood for CTCRW models
#' 
#' This function is designed for primary use within the \code{\link{crwMLE}}
#' model fitting function. But, it can be accessed for advanced \code{R} and
#' \code{crawl} users. Uses the state-space parameterization and Kalman filter
#' method presented in Johnson et al. (2008).
#' 

#' 
#' This function calls compiled C++ code which can be viewed in the
#' \code{src} directory of the crawl source package.
#' 
#' @param theta parameter values.
#' @param fixPar values of parameters held fixed (contains \code{NA} for
#' \code{theta} values).
#' @param y N by 2 matrix of coordinates with the longitude coordinate in the first column.
#' @param noObs vector with 1 for unobserved locations, and 0 for observed locations.
#' @param delta time difference to next location.
#' @param a initial state mean.
#' @param P intial state covariance matrix
#' @param mov.mf Movement covariate data.
#' @param err.mfX longitude error covariate data.
#' @param err.mfY latitude error covariate data.
#' @param rho A vector of known correlation coefficients for the error model, typically used for modern ARGOS data.
#' @param activity Stopping covariate (= 0 if animal is not moving).
#' @param n.errX number or longitude error parameters.
#' @param n.errY number of latitude error parameters.
#' @param n.mov number or movement parameters.
#' @param driftMod Logical. inicates whether a drift model is specified.
#' @param prior Function of theta that returns the log-density of the prior
#' @param need.hess Whether or not the Hessian will need to be calculated from
#' this call
#' @param constr Named list giving the parameter constraints
#' @return -2 * log-likelihood value for specified CTCRW model.
#' @author Devin S. Johnson
#' @seealso \code{\link{crwMLE}}
#' @references Johnson, D., J. London, M. -A. Lea, and J. Durban. 2008.
#' Continuous-time model for animal telemetry data. Ecology 89:1208-1215.
#' @export

crwN2ll = function(theta, fixPar, y, noObs, delta, a,
                   P, mov.mf, err.mfX, err.mfY, rho=NULL, activity=NULL,
                   n.errX, n.errY, n.mov, driftMod, prior, need.hess, 
                   constr=list(lower=-Inf, upper=Inf))
{
  if(!need.hess & any(theta < constr$lower | theta > constr$upper)) return(Inf)
  N <- nrow(y)
  par <- fixPar
  par[is.na(fixPar)] <- theta
  ###
  ### Process parameters for Fortran
  ###
  if (!is.null(err.mfX)) {
    theta.errX <- par[1:n.errX]
    Hmat <- exp(2 * err.mfX %*% theta.errX)
  } else Hmat <- rep(0.0, N)
  if (!is.null(err.mfY)) {
    theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
    Hmat <- cbind(Hmat,exp(2 * err.mfY %*% theta.errY))
  } else Hmat <- cbind(Hmat, Hmat)
  if(!is.null(rho)){
    Hmat = cbind(Hmat, sqrt(Hmat[,1])*sqrt(Hmat[,2])*rho)
  } else {Hmat = cbind(Hmat, rep(0,N))}
  Hmat[noObs==1,] = 0
  theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
  sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
  b <- exp(mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)])
  if (!is.null(activity)) {
    theta.stop <- par[(n.errX + n.errY + 2 * n.mov + 1)]
    b <- b / ((activity) ^ exp(theta.stop))
    active <- ifelse(b==Inf, 0, 1)
    b <- ifelse(b==Inf, 0, b) 
  } else {active=rep(1,N)}
  if (driftMod) {
    theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                         (n.errX + n.errY + 2 * n.mov + 2)]
    b.drift <- exp(log(b) - log(1+exp(theta.drift[2])))
    sig2.drift <- exp(log(sig2) + 2 * theta.drift[1]) 
    ll <- CTCRWNLL_DRIFT( y=as.matrix(y), Hmat, b, b.drift, sig2, sig2.drift, delta, noObs, active, a,  P)$ll
  } else {
    ll <- CTCRWNLL( y=as.matrix(y), Hmat, b, sig2, delta, noObs, active, a,  P)$ll
  }
  #movMats <- getQT(sig2, b, sig2.drift, b.drift, delta, driftMod=FALSE)

  if(is.null(prior)) return(-2 * ll)
  else return(-2 * (ll + prior(theta)))
}
