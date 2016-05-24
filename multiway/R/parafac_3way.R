parafac_3way <-
  function(data,nfac,xcx=sumsq(data),const=rep(0L,3),
           maxit=500,ctol=10^-4,Bfixed=NULL,Cfixed=NULL,
           Bstart=NULL,Cstart=NULL,Bstruc=NULL,Cstruc=NULL){
    # 3-way Parallel Factor Analysis (Parafac)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: January 27, 2016
    
    ### initialize Khatri-Rao product matrices
    xdims <- dim(data)
    if(is.null(Cfixed)){BkrA <- matrix(0,xdims[1]*xdims[2],nfac)}
    if(is.null(Bfixed)){CkrA <- matrix(0,xdims[1]*xdims[3],nfac)}
    CkrB <- matrix(0,xdims[2]*xdims[3],nfac)
    
    ### initialize reshaped data matrices
    Xa <- matrix(data,xdims[1],xdims[2]*xdims[3])
    if(is.null(Bfixed)){Xb <- matrix(aperm(data,perm=c(2,1,3)),xdims[2],xdims[1]*xdims[3])}
    if(is.null(Cfixed)){Xc <- matrix(aperm(data,perm=c(3,1,2)),xdims[3],xdims[1]*xdims[2])}
    rm(data)
    
    ### initialize parameter matrices
    if(const[1]==0L){
      Aold <- matrix(rnorm(xdims[1]*nfac),xdims[1],nfac)
    } else if(const[1]==1L){
      Aold <- svd(matrix(rnorm(xdims[1]*nfac),xdims[1],nfac),nu=nfac,nv=0)$u
    } else if(const[1]==2L){
      Aold <- Anew <- matrix(runif(xdims[1]*nfac),xdims[1],nfac)
    }
    if(is.null(Bfixed)){
      if(!is.null(Bstart)){
        Bold <- Bstart
      } else if(const[2]==0L){
        Bold <- matrix(rnorm(xdims[2]*nfac),xdims[2],nfac)
      } else if(const[2]==1L){
        Bold <- svd(matrix(rnorm(xdims[2]*nfac),xdims[2],nfac),nu=nfac,nv=0)$u
      } else if(const[2]==2L){
        Bold <- Bnew <- matrix(runif(xdims[2]*nfac),xdims[2],nfac)
      }
      if(!is.null(Bstruc)) Bold <- Bnew <- Bold * Bstruc
    } else {Bold <- Bnew <- Bfixed}
    if(is.null(Cfixed)){
      if(!is.null(Cstart)){
        Cold <- Cstart
      } else if(const[3]==0L){
        Cold <- matrix(rnorm(xdims[3]*nfac),xdims[3],nfac)
      } else if(const[3]==1L){
        Cold <- svd(matrix(rnorm(xdims[3]*nfac),xdims[3],nfac),nu=nfac,nv=0)$u
      } else if(const[3]==2L){
        Cold <- Cnew <- matrix(runif(xdims[3]*nfac),xdims[3],nfac)
      }
      if(!is.null(Cstruc)) Cold <- Cnew <- Cold * Cstruc
    } else {Cold <- Cnew <- Cfixed}
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      # Step 1: update mode A weights
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
      if(const[1]==0L){
        #Anew <- Xa%*%CkrB%*%smpower(crossprod(CkrB),-1)
        Anew <- Xa%*%CkrB%*%smpower(crossprod(Cold)*crossprod(Bold),-1)
      } else if(const[1]==1L) {
        Zmat <- Xa%*%CkrB
        Anew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
      } else if(const[1]==2L) {
        cpmat <- crossprod(CkrB)
        for(ii in 1:xdims[1]){Anew[ii,] <- fnnls(cpmat,crossprod(CkrB,Xa[ii,]))}
        if(any(colSums(Anew)==0)){
          Anew <- Aold
          vtol <- 0
          cflag <- 2
        }
      }
      
      # Step 2: update mode B weights
      if(is.null(Bfixed)){
        for(u in 1:nfac){CkrA[,u] <- kronecker(Cold[,u],Anew[,u])}
        if(is.null(Bstruc)){
          if(const[2]==0L){
            #Bnew <- Xb%*%CkrA%*%smpower(crossprod(CkrA),-1)
            Bnew <- Xb%*%CkrA%*%smpower(crossprod(Cold)*crossprod(Anew),-1)
          } else if(const[2]==1L) {
            Zmat <- Xb%*%CkrA
            Bnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[2]==2L) {
            cpmat <- crossprod(CkrA)
            for(ii in 1:xdims[2]){Bnew[ii,] <- fnnls(cpmat,crossprod(CkrA,Xb[ii,]))}
            if(any(colSums(Bnew)==0)){
              Bnew <- Bold
              vtol <- 0
              cflag <- 2
            }
          } # end if(const[2]==0L)
        } else {
          for(u in 1:nfac){
            Zhat = Xb - tcrossprod(Bnew[,-u],CkrA[,-u])
            Bnew[,u] = ( (Zhat %*% CkrA[,u]) / sum(CkrA[,u]^2) ) * Bstruc[,u]
          }
        } # end if(is.null(Bstruc))
      } # end if(is.null(Bfixed))
      
      # Step 3: update mode C weights
      if(is.null(Cfixed)){
        for(u in 1:nfac){BkrA[,u] <- kronecker(Bnew[,u],Anew[,u])}
        if(is.null(Cstruc)){
          if(const[3]==0L){
            #Cnew <- Xc%*%BkrA%*%smpower(crossprod(BkrA),-1)
            Cnew <- Xc%*%BkrA%*%smpower(crossprod(Bnew)*crossprod(Anew),-1)
          } else if(const[3]==1L) {
            Zmat <- Xc%*%BkrA
            Cnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[3]==2L) {
            cpmat <- crossprod(BkrA)
            for(ii in 1:xdims[3]){Cnew[ii,] <- fnnls(cpmat,crossprod(BkrA,Xc[ii,]))}
            if(any(colSums(Cnew)==0)){
              Cnew <- Cold
              vtol <- 0
              cflag <- 2
            }
          } # end if(const[3]==0L)
        } else {
          for(u in 1:nfac){
            Zhat = Xc - tcrossprod(Cnew[,-u],BkrA[,-u])
            Cnew[,u] = ( (Zhat %*% BkrA[,u]) / sum(BkrA[,u]^2) ) * Cstruc[,u]
          }
        } # end if(is.null(Cstruc))
      } # end if(is.null(Cfixed))
      
      # Step 4: check for convergence
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cnew[,u],Bnew[,u])}
      ssenew <- sum((Xa-tcrossprod(Anew,CkrB))^2)
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Aold <- Anew
      Bold <- Bnew
      Cold <- Cnew
      sseold <- ssenew
      iter <- iter + 1
      
    } # end while(vtol>ctol && iter<maxit)
    
    ### scale and order solution
    if(is.null(Bfixed) & is.null(Cfixed)){
      
      # put the scale in Mode C
      adg <- colMeans(Anew^2)
      Anew <- Anew%*%(diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bnew^2)
      Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
      Cnew <- Cnew%*%(diag(nfac)*((adg*bdg)^0.5))
      
      # order according to sum-of-squares
      fordr <- order(colSums(Cnew^2),decreasing=TRUE)
      Anew <- as.matrix(Anew[,fordr])
      Bnew <- as.matrix(Bnew[,fordr])
      Cnew <- as.matrix(Cnew[,fordr])
      
    }
    
    ### GCV criterion
    Adf <- nfac * ifelse(const[1]==1L, xdims[1]-(nfac+1)/2, xdims[1]-1L)
    Bdf <- ifelse(is.null(Bstruc),
                  nfac * ifelse(const[2]==1L, xdims[2]-(nfac+1)/2, xdims[2]-1L),
                  sum(Bstruc) - nfac)
    Cdf <- ifelse(is.null(Cstruc),
                  nfac * ifelse(const[3]==1L, xdims[3]-(nfac-1)/2, xdims[3]),
                  sum(Cstruc))
    edf <- c(Adf,Bdf,Cdf)
    pxdim <- prod(xdims)
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C")
    pfac <- list(A=Anew,B=Bnew,C=Cnew,Rsq=Rsq,GCV=GCV,edf=edf,
                 iter=iter,cflag=cflag,const=const)
    return(pfac)
    
  }