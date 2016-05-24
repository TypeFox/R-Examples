parafac2_3way <- 
  function(data,nfac,xcx=sumsq(data),const=rep(0L,3),
           maxit=500,ctol=10^-4,Gfixed=NULL,Bfixed=NULL,
           Cfixed=NULL,Gstart=NULL,Bstart=NULL,Cstart=NULL,
           Gstruc=NULL,Bstruc=NULL,Cstruc=NULL){
    # 3-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 15, 2016
    
    ### initialize Khatri-Rao product matrices
    xdims <- rep(NA,3)
    xdims[2] <- ncol(data[[1]])
    xdims[3] <- length(data)
    if(is.null(Cfixed)){BkrA <- matrix(0,nfac*xdims[2],nfac)}
    if(is.null(Bfixed)){CkrA <- matrix(0,nfac*xdims[3],nfac)}
    CkrB <- matrix(0,xdims[2]*xdims[3],nfac)
    
    ### initialize stuff for Mode A update
    Rknew <- vector("list",xdims[3])
    Xtilde <- array(0,dim=c(nfac,xdims[2],xdims[3]))
    
    ### initialize parameter matrices
    if(is.null(Gfixed)){
      if(!is.null(Gstart)){
        Gold <- Gstart
      } else if(const[1]==0L){
        Gold <- matrix(rnorm(nfac^2),nfac,nfac)
        #Gold <- Gold %*% (diag(nfac)/sqrt(colSums(Gold^2)))
      } else if(const[1]==1L){
        Gold <- Gnew <- diag(nfac)
      }
      if(!is.null(Gstruc)) Gold <- Gnew <- Gold * Gstruc
    } else {Gold <- Gnew <- Gfixed}
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
        #Cold <- matrix(rnorm(xdims[3]*nfac),xdims[3],nfac)
        Cold <- matrix(runif(xdims[3]*nfac),xdims[3],nfac)
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
      
      ## Step 1: update mode A weights
      # 1a: update orthogonal projections
      for(kk in 1:xdims[3]){
        xsvd <- svd(data[[kk]]%*%Bold%*%tcrossprod((diag(nfac)*Cold[kk,]),Gold))
        Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
        Xtilde[,,kk] <- crossprod(Rknew[[kk]],data[[kk]])
      }
      # 1b: update correlation matrix
      if(is.null(Gfixed)){
        if(const[1]==0L){
          Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3])
          for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
          if(is.null(Gstruc)){
            #Gnew <- Xa%*%CkrB%*%smpower(crossprod(CkrB),-1)
            Gnew <- Xa%*%CkrB%*%smpower(crossprod(Cold)*crossprod(Bold),-1)
          } else {
            for(u in 1:nfac){
              Zhat = Xa - tcrossprod(Gnew[,-u],CkrB[,-u])
              Gnew[,u] = ( (Zhat %*% CkrB[,u]) / sum(CkrB[,u]^2) ) * Gstruc[,u]
            }
          } # end if(is.null(Gstruc))
        } # end if(const[1]==0L)
      } # end if(is.null(Gfixed))
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed)){
        Xb <- matrix(aperm(Xtilde,perm=c(2,1,3)),xdims[2],nfac*xdims[3])
        for(u in 1:nfac){CkrA[,u] <- kronecker(Cold[,u],Gnew[,u])}
        if(is.null(Bstruc)){
          if(const[2]==0L){
            #Bnew <- Xb%*%CkrA%*%smpower(crossprod(CkrA),-1)
            Bnew <- Xb%*%CkrA%*%smpower(crossprod(Cold)*crossprod(Gnew),-1)
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
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed)){
        Xc <- matrix(aperm(Xtilde,perm=c(3,1,2)),xdims[3],nfac*xdims[2])
        for(u in 1:nfac){BkrA[,u] <- kronecker(Bnew[,u],Gnew[,u])}
        if(is.null(Cstruc)){
          if(const[3]==0L){
            #Cnew <- Xc%*%BkrA%*%smpower(crossprod(BkrA),-1)
            Cnew <- Xc%*%BkrA%*%smpower(crossprod(Bnew)*crossprod(Gnew),-1)
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
      
      ## Step 4: check for convergence
      ssenew <- 0
      for(kk in 1:xdims[3]){
        ssenew <- ssenew + sum((data[[kk]]-tcrossprod(Rknew[[kk]]%*%Gnew%*%(diag(nfac)*Cnew[kk,]),Bnew))^2)
      }
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Gold <- Gnew
      Bold <- Bnew
      Cold <- Cnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### scale and order solution
    if(is.null(Gfixed) & is.null(Bfixed) & is.null(Cfixed)){
      
      # put the scale in Mode C
      adg <- colSums(Gnew^2)
      Gnew <- Gnew%*%(diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bnew^2)
      Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
      Cnew <- Cnew%*%(diag(nfac)*((adg*bdg)^0.5))
      
      # order according to sum-of-squares
      fordr <- order(colSums(Cnew^2),decreasing=TRUE)
      Gnew <- as.matrix(Gnew[,fordr])
      Bnew <- as.matrix(Bnew[,fordr])
      Cnew <- as.matrix(Cnew[,fordr])
      
    }
    
    ### GCV criterion
    ntotal <- sum(sapply(data,nrow))
    Adf <- nfac * (ntotal - xdims[3]*(nfac+1)/2)
    if(!is.null(Gstruc)) GtG <- crossprod(Gstruc)
    Gdf <- ifelse(is.null(Gstruc),
                  ifelse(const[1]==1L, 0, nfac*(nfac-1)/2),
                  sum(GtG[lower.tri(GtG)]>0L))
    Bdf <- ifelse(is.null(Bstruc),
                  nfac * ifelse(const[2]==1L, xdims[2]-(nfac+1)/2, xdims[2]-1L),
                  sum(Bstruc) - nfac)
    Cdf <- ifelse(is.null(Cstruc),
                  nfac * ifelse(const[3]==1L, xdims[3]-(nfac-1)/2, xdims[3]),
                  sum(Cstruc))
    edf <- c(Adf+Gdf,Bdf,Cdf)
    pxdim <- ntotal * xdims[2]
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C")
    pfac <- list(A=list(H=Rknew,G=Gnew),B=Bnew,C=Cnew,Rsq=Rsq,GCV=GCV,
                 edf=edf,iter=iter,cflag=cflag,const=const)
    return(pfac)
    
  }