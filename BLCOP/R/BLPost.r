###############################################################################
# Mango Solutions
# posteriorEst
# Author: Francisco
# $Rev: 4763 $
# $LastChangedDate: 2010-02-23 20:32:39 +0000 (Tue, 23 Feb 2010) $
#
###############################################################################
# DESCRIPTION: Computes the Black-Litterman posterior estimate 
# KEYWORDS: math
###############################################################################

#' This function performs the "core" calculation of the Black-Litterman model.  
#' @param views An object of class BLViews
#' @param mu A vector of mean equilibrium returns 
#' @param tau The "tau" parameter in the Black-Litterman model.
#' @param sigma The variance-covariance matrix of the returns of the assets
#' @param kappa if greater than 0, the confidences in each view are replaced.  See the online help for details
#' @return An object of class BLResult holding the updated Black-Litterman posterior
#' @author Francisco
#' @export

posteriorEst <- function
(
    views,     # full BLview object.  
	mu,    # Equilibrium expected returns
	tau = 0.5,       # Degree of uncertainty in prior
    
    sigma,     # variance-covariance matrix of asset returns
    kappa = 0  # if greater than 0, the view confidences will be ignored and the
               # omega matrix in the BL model will be replaced by kappa * P %*% sigma %*% t(P)
)
{

  # preallocate the return
  numAssets <- length(assetSet(views))
  
  P <- views@P
  if(kappa == 0)
  {    
      if(length(views@confidences) > 1) 
        omega <- diag( 1/ views@confidences) 
      else 
        omega <- matrix(1 / views@confidences, 1,1)
  }
  
  else    
  {
      omega <- kappa * tcrossprod(P %*% sigma, P)
      omegaInv <- solve(omega)
  }     

  qv <- views@qv
  sigmaInv <- solve(sigma)

  # The following steps are the core Black-Litterman calculations
  
  temp <- tcrossprod(sigma, P)
 
  postMu <- mu + tau * temp %*% solve(tau * P %*% temp + omega, qv - P %*% mu)
  postMu <- as.numeric(postMu)
  
  postSigma <- (1 + tau) * sigma - tau^2 * temp %*% solve(tau * P %*% temp + omega, P %*% sigma)
  names(mu) <- assetSet(views)
  names(postMu) <- assetSet(views)
  
  new("BLResult", views = views, tau = tau, priorMean = mu, priorCovar = sigma,
               posteriorMean = postMu, posteriorCovar = postSigma, kappa = kappa )
}

#' BLposterior
#' @param returns A matrix of time series of returns.  The columns should correspond to individual assets.
#' @param views An object of class BLViews
#' @param tau The "tau" parameter in the Black-Litterman model.
#' @param marketIndex A set of returns of a market index.
#' @param riskFree A time series of risk-free rates of return.  Defaults to 0
#' @param kappa if greater than 0, the confidences in each view are replaced.  See the online help for details  
#' @param covEstimator A string holding the name of the function that should be used to estimate the variance-covariance matrix.
#'     This function should simply return a matrix.
#' @return An object of class BLResult
#' @author Francisco
#' @export

BLPosterior <- function
(
  returns,
  views,
  tau = 1,         
  marketIndex,
  riskFree = NULL,
  kappa = 0,
  covEstimator = "cov"
)
{
  covEstimator <- match.fun(covEstimator)
  alphaInfo <- CAPMList(returns, marketIndex, riskFree = riskFree)
  post <- posteriorEst(views, tau = tau, mu = alphaInfo[["alphas"]], 
      sigma = unclass(covEstimator(returns)),  kappa = kappa)
  post
}