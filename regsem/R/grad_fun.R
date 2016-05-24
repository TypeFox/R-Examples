
grad_fun = function(par,ImpCov,SampCov,Areg,Sreg,A,S,F,
                    A_fixed,A_est,S_fixed,S_est,lambda,type,pars_pen,I){


  # m = dim(ImpCov)[1]
  #grad = rep(0,length(par))
  #facLoad <- matrix(c(1,par[1:(nload-1)]),ncol(A)-1,nfac) #extract just the factor loadings
  #facLoad <- A[1:nload,ncol(A)]
  #facLoad <- A[1:nvar,(ncol(A)-nfac +1):ncol(A)]
  #facCov <- S[(nrow(S)-nfac + 1):nrow(S),(ncol(S)-nfac + 1):ncol(S)]
  #onemat <- matrix(1,nrow(facLoad),ncol(facLoad)) # matrix of ones for fac load deriv

  # from Cudeck 1993
  #h <- 0_00001
  # all possible parameters
  #pPos <- nvar*nfac + nfac*(nfac + 1)/2 + nvar*(nvar + 1)/2
  #pPos <- rep(1,length(par))

   grad_out <- rep(0,length(par))
  # h <- 0.00001
   h = sqrt(.Machine$double.eps)

   A_iter <- max(A)


  if(type=="none"){

    for(i in 1:length(par)){
      add <- rep(0,length(par))
      add[i] <- h
      ImpCovL = rcpp_RAMmult((par+add),A,S,S_fixed,A_fixed,A_est,S_est,F,I)[[1]]
      #ImpCovL = RAMmult((par+add),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
      ImpCovDot <- (ImpCovL - ImpCov)/h
      grad_out[i] <- 0.5 * trace(solve(ImpCov) %*% (ImpCov - SampCov) %*% solve(ImpCov) %*% ImpCovDot)
    }


  }  else if(type=="ridge"){ # not specified

    for(i in 1:length(par)){
      add <- rep(0,length(par))
      add[i] <- h
      ImpCovL = rcpp_RAMmult((par+add),A,S,S_fixed,A_fixed,A_est,S_est,F,I)[[1]]
      #ImpCovL = RAMmult((par+add),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
      ImpCovDot <- (ImpCovL - ImpCov)/h
      grad_out[i] <- 0.5 * (trace(solve(ImpCov) %*% (ImpCov - SampCov) %*% solve(ImpCov) %*% ImpCovDot)) +
        if(any(i==pars_pen)) 2*lambda*(max(Areg[A == i], Sreg[S==i])) else(0)
    }


  }  else if(type=="lasso"){
    #lasso

    for(i in 1:length(par)){
      add <- rep(0,length(par))
      add[i] <- h
      ImpCovL = rcpp_RAMmult((par+add),A,S,S_fixed,A_fixed,A_est,S_est,F,I)[[1]]
      #ImpCovL = RAMmult((par+add),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
      ImpCovDot <- (ImpCovL - ImpCov)/h
      grad_out[i] <- 0.5 * (trace(solve(ImpCov) %*% (ImpCov - SampCov) %*% solve(ImpCov) %*% ImpCovDot)) +
        if(any(i==pars_pen)) lambda*sign(max(Areg[A == i], Sreg[S==i])) else(0)
    }


  }  else if(type=="enet"){ # not specified
    #elastic net

    for(i in 1:length(par)){
      add <- rep(0,length(par))
      add[i] <- 1 + h
      ImpCovL = RAMmult((par + add),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
      ImpCovDot <- (ImpCovL - ImpCov)/h
      grad_out[i] <- 0.5 * trace(solve(ImpCov) %*% (ImpCov - SampCov) %*% solve(ImpCov) %*% ImpCovDot)
    }


  }else if(type=="ols_lasso"){
    #lasso

    for(i in 1:length(par)){
      add <- rep(0,length(par))
      add[i] <- h
      ImpCovL = RAMmult((par + add),A,S,F,A_fixed,A_est,S_fixed,S_est)[[1]]
      ImpCovDot <- (ImpCovL - ImpCov)/h
      grad_out[i] <- 0.5 * trace(solve(ImpCov) %*% (ImpCov - SampCov) %*%
                                   solve(ImpCov) %*% ImpCovDot) + if(i <= A_iter) lambda*sign(Areg[A==i]) else(0)
    }


  }

  as.numeric(grad_out)


}
