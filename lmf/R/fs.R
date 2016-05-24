fs <-
function(At,
               at,
               npar,
               nyear,
               method,
               control,
               ...)
{
  #Estimate average alpha covariance matrix components, calculate
  #the upper triangular matrix of the cholesky decomposition and extract
  #the components to a vector
  D <- as.vector(chol((Reduce('+', At) / nyear)))
  #Remove the zero elements
  D <- D[D != 0]
  #Optimize the log likelihood function with respect to M
  D.opt <- optim(par = D, fn = lnL.M, gr = NULL, At = At,
                 at = at, npar = npar, method = method,
                 control = control, ...)
  #Create a matrix of equal size as the M matrix
  Dm <- matrix(rep(0, npar^2), ncol = npar)
  #Insert the optimized elements of D into Dm to create the cholesky
  #decomposed optimal upper triangular matrix
  Dm[upper.tri(Dm, diag = TRUE)] <- D.opt$par
  #Replace any negative diagonal elements of Dm by their absolute value
  diag(Dm) <- abs(diag(Dm))
  #Extract the convergence indicator from optim and assigne "yes" or "no"
  ret <- list(convergence = ifelse(D.opt$convergence == 0, "yes", "no"))
  #Extract the number of iterations from optim
  ret$iterations <- as.numeric(D.opt$counts)[1] * (npar * (npar + 1) / 2)
  #Calculate the temporal covariance matrix, M, if function converged
  if(ret$convergence == "yes")
  {
    ret$M <- crossprod(Dm)
    #If all eigen values is not positive (or is complex) the estimated matrix
    #is not positive definite and a "close" positive definite matrix
    #is returned in its place through the algorithm in the function nearPD
    if(is.complex(eigen(ret$M)$values > 0))
      ret$M <- nearPD(ret$M)
    if(!all(eigen(ret$M)$values > 0))
      ret$M <- nearPD(ret$M)
    #Estimate the alphas with the optimal covariance matrix, M
    ret$aM <- lnL.M(D = D.opt$par, At = At, at = at, npar = npar,
                    ret.alphas = TRUE)
  }
  else{
    # Return empty matrix (filled with NA)
    ret$M <- matrix(rep(NA, npar^2), ncol = npar)
    # Return empty vector (filled with NA)
    ret$aM <- rep(NA, npar)
    # Also return warning
    warning("optim() did not reach convergence. Temporal covariance matrix (M) and temporal mean coefficients a(M) unavailable")
  }
  #Add row and column names
  dimnames(ret$M) <- dimnames(At[[1]])
  #Output
  ret
}
