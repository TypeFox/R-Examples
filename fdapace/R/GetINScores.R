# This function obtains the FPC scores for dense
# regular functional data by trapezoidal rule integration

######
# Input:
######  
# ymat: n by p matrix of dense regular functional observations 
# t: list of observed time grids for the functional observations
######
# Output: 
######
# ret: a list of:
#        xiEst: n by length(lambda) matrix of estimated FPC scores
#        fittedY: n by p matrix of fitted/recovered functional observations
##########################################################################

GetINScores <- function(ymat, t, optns, mu, lambda, phi, sigma2=NULL){
  if(length(lambda) != ncol(phi)){
    stop('No. of eigenvalues is not the same as the no. of eigenfunctions.')
  }

  n = nrow(ymat)
  tau = sort(unique(signif( unlist(t),14 ))) # get observed time grid
  ranget <- diff(range(tau))
  mumat = matrix(rep(mu, n), nrow = n, byrow = TRUE)
  cymat = ymat - mumat

  xiEst = matrix(0, nrow = n, ncol = length(lambda)) 
  # Get Scores xiEst
  for(i in 1:length(lambda)){
    tempmat = cymat * matrix(rep(phi[,i],n), nrow = n, byrow = TRUE)
    xiEst[,i] = sapply(1:n, function(j) trapzRcpp(X = tau[!is.na(tempmat[j,])], Y = tempmat[j, !is.na(tempmat[j,])]))
    if (optns[['shrink']] && !is.null(sigma2)) {
      xiEst[, i] <- xiEst[, i] * lambda[i] / 
                    (lambda[i] + ranget * sigma2 / length(tau))
    }
  }

  # Get Fitted Y: n by p matrix on observed time grid
  fittedY = mumat + t(phi %*% t(xiEst))

  ret = list('xiEst' = xiEst, xiVar = NULL, 'fittedY' = fittedY)

  return(ret)
}
