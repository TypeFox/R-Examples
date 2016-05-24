
hessian = function(par,ImpCov,SampCov,A,A_fixed,A_est,S,S_fixed,S_est,F,lambda,alpha,type){

  hess_out <- matrix(0,length(par),length(par))
  h = sqrt(.Machine$double.eps)
  #h = .000001

  # using Cudeck, Klebe, Henly (1993)
  # p_113 Eq_ 1 (Magnus & Neudecker, 2007)

  add <- matrix(0,length(par),length(par))
  diag(add) <- h

  for(i in 1:nrow(hess_out)){
    for(j in 1:ncol(hess_out)){
      ImpCovI = RAMmult((par + add[i,]),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
      ImpCovII <- (ImpCovI - ImpCov)/h
      ImpCovJ = RAMmult((par + add[,j]),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
      ImpCovJJ <- (ImpCovJ - ImpCov)/h
      hess_out[i,j] <- 0.5 * trace(solve(ImpCov) %*% ImpCovII %*% solve(ImpCov) %*% ImpCovJJ)
    }
  }



  hess_out


}

