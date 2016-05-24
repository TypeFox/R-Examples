tucker_3way <-
  function(data,nfac,xcx=sumsq(data),maxit=500,
           ctol=10^-4,Afixed=NULL,Bfixed=NULL,Cfixed=NULL,
           Bstart=NULL,Cstart=NULL){
    # 3-way Tucker model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: October 14, 2015
    
    ### initialize reshaped data matrices
    xdims <- dim(data)
    Xa <- matrix(data,xdims[1],xdims[2]*xdims[3])
    if(is.null(Bfixed)){Xb <- matrix(aperm(data,perm=c(2,1,3)),xdims[2],xdims[1]*xdims[3])}
    if(is.null(Cfixed)){Xc <- matrix(aperm(data,perm=c(3,1,2)),xdims[3],xdims[1]*xdims[2])}
    rm(data)
    
    ### initialize parameter matrices
    if(is.null(Afixed)){
      Aold <- svd(matrix(rnorm(xdims[1]*nfac[1]),xdims[1],nfac[1]),nu=nfac[1],nv=0)$u
    } else {Aold <- Anew <- Afixed}
    if(is.null(Bfixed)){
      if(is.null(Bstart)){
        Bold <- svd(matrix(rnorm(xdims[2]*nfac[2]),xdims[2],nfac[2]),nu=nfac[2],nv=0)$u
      } else {Bold <- Bstart}
    } else {Bold <- Bnew <- Bfixed}
    if(is.null(Cfixed)){
      if(is.null(Cstart)){
        Cold <- svd(matrix(rnorm(xdims[3]*nfac[3]),xdims[3],nfac[3]),nu=nfac[3],nv=0)$u
      } else {Cold <- Cstart}
    } else {Cold <- Cnew <- Cfixed}
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      # Step 1: update mode A weights
      if(is.null(Afixed)){ Anew <- svd(Xa%*%kronecker(Cold,Bold),nu=nfac[1],nv=0)$u }
      
      # Step 2: update mode B weights
      if(is.null(Bfixed)){ Bnew <- svd(Xb%*%kronecker(Cold,Anew),nu=nfac[2],nv=0)$u }
      
      # Step 3: update mode C weights
      if(is.null(Cfixed)){ Cnew <- svd(Xc%*%kronecker(Bnew,Anew),nu=nfac[3],nv=0)$u }
      
      # Step 4: update G and check for convergence
      Ga <- crossprod(Anew,Xa%*%kronecker(Cnew,Bnew))
      ssenew <- xcx - sum(rowSums(Ga^2))
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Aold <- Anew
      Bold <- Bnew
      Cold <- Cnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### GCV criterion
    Adf <- xdims[1]*nfac[1] - nfac[1]*(nfac[1]+1)/2
    Bdf <- xdims[2]*nfac[2] - nfac[2]*(nfac[2]+1)/2
    Cdf <- xdims[3]*nfac[3] - nfac[3]*(nfac[3]+1)/2
    edf <- c(Adf,Bdf,Cdf,prod(nfac))
    pxdim <- prod(xdims)
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C","G")
    tuck <- list(A=Anew,B=Bnew,C=Cnew,G=array(Ga,dim=nfac),
                 Rsq=Rsq,GCV=GCV,edf=edf,iter=iter,cflag=cflag)
    return(tuck)
    
  }