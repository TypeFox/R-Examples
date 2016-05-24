# Adapting the Cholesky decomposition of the covariance matrix
# using the approach of Vihola (2011)
# 
# @param nPar Number of parameters.
# @param iter Index of the current MCMC iteration.
# @param S    Matrix containig the transpose of the Cholesky decomposition of the 
#             covariance matrix of the proposal.
# @param U    Vector of normal perturbations used to perturb the parameters (i.e. rnorm(nPar)).
# @param gamma Controls the speed of adaption. Should be between 0.5 and 1. 
#              A lower gamma leads to faster adaption.
# @param alpha Acceptance ratio. Typically min(1, propLikelihood / oldLikelihood).
# @param accRate Targe acceptance rate
# @return An matrix containing adapted version of transpose of the 
#         Cholesky decomposition of the covariance matrix of the proposal.
# @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> adapted from Vihola's adaptMCMC package                       
# @export Not exporting it at the moment
# 
.adaptChol <- function(nPar, iter, S, U, gamma, alpha, accRate) 
{  
  adaptRate <- min(5, nPar * iter^(-gamma))
  M <- S %*% (diag(nPar) + adaptRate * (alpha - accRate) * U %*% t(U)/sum(U^2)) %*% t(S)
  eig <- eigen(M, only.values = TRUE)$values
  tol <- ncol(M) * max(abs(eig)) * .Machine$double.eps
  
  if (!isSymmetric(M) | is.complex(eig) | !all(Re(eig) >  tol)) 
  { 
    M <- as.matrix(nearPD(M)$mat)
  }
  
  # Rounding at 10e-14 otherwise there is a feedback mechanism 
  # and the results are not reproducible
  S <- round( t(chol(M)), 14)
  
  return(S)
}
