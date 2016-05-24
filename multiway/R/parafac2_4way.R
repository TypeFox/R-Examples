parafac2_4way <- 
  function(data,nfac,xcx=sumsq(data),const=rep(0L,4),
           maxit=500,ctol=10^-4,Gfixed=NULL,Bfixed=NULL,Cfixed=NULL,
           Dfixed=NULL,Gstart=NULL,Bstart=NULL,Cstart=NULL,Dstart=NULL,
           Gstruc=NULL,Bstruc=NULL,Cstruc=NULL,Dstruc=NULL){
    # 4-way Parallel Factor Analysis 2 (Parafac2)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 15, 2015
    
    ### initialize Khatri-Rao product matrices
    xdims <- rep(NA,4)
    xdims[2] <- dim(data[[1]])[2]
    xdims[3] <- dim(data[[1]])[3]
    xdims[4] <- length(data)
    if(is.null(Dfixed)){CBkrA <- matrix(0,nfac*xdims[2]*xdims[3],nfac)}
    if(is.null(Cfixed)){DBkrA <- matrix(0,nfac*xdims[2]*xdims[4],nfac)}
    if(is.null(Bfixed)){DCkrA <- matrix(0,nfac*xdims[3]*xdims[4],nfac)}
    DCkrB <- matrix(0,xdims[2]*xdims[3]*xdims[4],nfac)
    CkrB <- matrix(0,xdims[2]*xdims[3],nfac)
    
    ### initialize stuff for Mode A update
    Rknew <- vector("list",xdims[4])
    Xtilde <- array(0,dim=c(nfac,xdims[2],xdims[3],xdims[4]))
    
    ### reshape raw data
    for(kk in 1:xdims[4]){
      mdim <- dim(data[[kk]])
      data[[kk]] <- matrix(data[[kk]],mdim[1],xdims[2]*xdims[3])
    }
    
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
        Cold <- matrix(rnorm(xdims[3]*nfac),xdims[3],nfac)
      } else if(const[3]==1L){
        Cold <- svd(matrix(rnorm(xdims[3]*nfac),xdims[3],nfac),nu=nfac,nv=0)$u
      } else if(const[3]==2L){
        Cold <- Cnew <- matrix(runif(xdims[3]*nfac),xdims[3],nfac)
      }
      if(!is.null(Cstruc)) Cold <- Cnew <- Cold * Cstruc
    } else {Cold <- Cnew <- Cfixed}
    if(is.null(Dfixed)){
      if(!is.null(Dstart)){
        Dold <- Dstart
      } else if(const[4]==0L){
        #Dold <- matrix(rnorm(xdims[4]*nfac),xdims[4],nfac)
        Dold <- matrix(runif(xdims[4]*nfac),xdims[4],nfac)
      } else if(const[4]==1L){
        Dold <- svd(matrix(rnorm(xdims[4]*nfac),xdims[4],nfac),nu=nfac,nv=0)$u
      } else if(const[4]==2L){
        Dold <- Dnew <- matrix(runif(xdims[4]*nfac),xdims[4],nfac)
      }
      if(!is.null(Dstruc)) Dold <- Dnew <- Dold * Dstruc
    } else {Dold <- Dnew <- Dfixed}
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      ## Step 1: update mode A weights
      # 1a: update orthogonal projections
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cold[,u],Bold[,u])}
      for(kk in 1:xdims[4]){
        xsvd <- svd(data[[kk]]%*%CkrB%*%tcrossprod((diag(nfac)*Dold[kk,]),Gold))
        Rknew[[kk]] <- tcrossprod(xsvd$u,xsvd$v)
        Xtilde[,,,kk] <- array(crossprod(Rknew[[kk]],data[[kk]]),dim=c(nfac,xdims[2],xdims[3]))
      }
      # 1b: update correlation matrix
      if(is.null(Gfixed)){
        if(const[1]==0L){
          Xa <- matrix(Xtilde,nfac,xdims[2]*xdims[3]*xdims[4])
          for(u in 1:nfac){DCkrB[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Bold[,u]))}
          if(is.null(Gstruc)){
            #Gnew <- Xa%*%DCkrB%*%smpower(crossprod(DCkrB),-1)
            Gnew <- Xa%*%DCkrB%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Bold),-1)
          } else {
            for(u in 1:nfac){
              Zhat = Xa - tcrossprod(Gnew[,-u],DCkrB[,-u])
              Gnew[,u] = ( (Zhat %*% DCkrB[,u]) / sum(DCkrB[,u]^2) ) * Gstruc[,u]
            }
          } # end if(is.null(Gstruc))
        } # end if(const[1]==0L)
      } # end if(is.null(Gfixed))
      
      ## Step 2: update mode B weights
      if(is.null(Bfixed)){
        Xb <- matrix(aperm(Xtilde,perm=c(2,1,3,4)),xdims[2],nfac*xdims[3]*xdims[4])
        for(u in 1:nfac){DCkrA[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Gnew[,u]))}
        if(is.null(Bstruc)){
          if(const[2]==0L){
            #Bnew <- Xb%*%DCkrA%*%smpower(crossprod(DCkrA),-1)
            Bnew <- Xb%*%DCkrA%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Gnew),-1)
          } else if(const[2]==1L) {
            Zmat <- Xb%*%DCkrA
            Bnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[2]==2L) {
            cpmat <- crossprod(DCkrA)
            for(ii in 1:xdims[2]){Bnew[ii,] <- fnnls(cpmat,crossprod(DCkrA,Xb[ii,]))}
            if(any(colSums(Bnew)==0)){
              Bnew <- Bold
              vtol <- 0
              cflag <- 2
            }
          } # end if(const[2]==0L)
        } else {
          for(u in 1:nfac){
            Zhat = Xb - tcrossprod(Bnew[,-u],DCkrA[,-u])
            Bnew[,u] = ( (Zhat %*% DCkrA[,u]) / sum(DCkrA[,u]^2) ) * Bstruc[,u]
          }
        } # end if(is.null(Bstruc))
      } # end if(is.null(Bfixed))
      
      ## Step 3: update mode C weights
      if(is.null(Cfixed)){
        Xc <- matrix(aperm(Xtilde,perm=c(3,1,2,4)),xdims[3],nfac*xdims[2]*xdims[4])
        for(u in 1:nfac){DBkrA[,u] <- kronecker(Dold[,u],kronecker(Bnew[,u],Gnew[,u]))}
        if(is.null(Cstruc)){
          if(const[3]==0L){
            #Cnew <- Xc%*%DBkrA%*%smpower(crossprod(DBkrA),-1)
            Cnew <- Xc%*%DBkrA%*%smpower(crossprod(Dold)*crossprod(Bnew)*crossprod(Gnew),-1)
          } else if(const[3]==1L) {
            Zmat <- Xc%*%DBkrA
            Cnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[3]==2L) {
            cpmat <- crossprod(DBkrA)
            for(ii in 1:xdims[3]){Cnew[ii,] <- fnnls(cpmat,crossprod(DBkrA,Xc[ii,]))}
            if(any(colSums(Cnew)==0)){
              Cnew <- Cold
              vtol <- 0
              cflag <- 2
            }
          } # end if(const[3]==0L)
        } else {
          for(u in 1:nfac){
            Zhat = Xc - tcrossprod(Cnew[,-u],DBkrA[,-u])
            Cnew[,u] = ( (Zhat %*% DBkrA[,u]) / sum(DBkrA[,u]^2) ) * Cstruc[,u]
          }
        } # end if(is.null(Cstruc))
      } # end if(is.null(Cfixed))
      
      ## Step 4: update mode D weights
      if(is.null(Dfixed)){
        Xd <- matrix(aperm(Xtilde,perm=c(4,1,2,3)),xdims[4],nfac*xdims[2]*xdims[3])
        for(u in 1:nfac){CBkrA[,u] <- kronecker(Cnew[,u],kronecker(Bnew[,u],Gnew[,u]))}
        if(is.null(Dstruc)){
          if(const[4]==0L){
            #Dnew <- Xd%*%CBkrA%*%smpower(crossprod(CBkrA),-1)
            Dnew <- Xd%*%CBkrA%*%smpower(crossprod(Cnew)*crossprod(Bnew)*crossprod(Gnew),-1)
          } else if(const[4]==1L) {
            Zmat <- Xd%*%CBkrA
            Dnew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
          } else if(const[4]==2L) {
            cpmat <- crossprod(CBkrA)
            for(ii in 1:xdims[4]){Dnew[ii,] <- fnnls(cpmat,crossprod(CBkrA,Xd[ii,]))}
            if(any(colSums(Dnew)==0)){
              Dnew <- Dold
              vtol <- 0
              cflag <- 2
            }
          } # end if(const[4]==0L)
        } else {
          for(u in 1:nfac){
            Zhat = Xd - tcrossprod(Dnew[,-u],CBkrA[,-u])
            Dnew[,u] = ( (Zhat %*% CBkrA[,u]) / sum(CBkrA[,u]^2) ) * Dstruc[,u]
          }
        } # end if(is.null(Dstruc))
      } # end if(is.null(Dfixed))
      
      ## Step 5: check for convergence
      for(u in 1:nfac){CkrB[,u] <- kronecker(Cnew[,u],Bnew[,u])}
      ssenew <- 0
      for(kk in 1:xdims[4]){
        ssenew <- ssenew + sum((data[[kk]]-tcrossprod(Rknew[[kk]]%*%Gnew%*%(diag(nfac)*Dnew[kk,]),CkrB))^2)
      }
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Gold <- Gnew
      Bold <- Bnew
      Cold <- Cnew
      Dold <- Dnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### scale and order solution
    if(is.null(Gfixed) & is.null(Bfixed) & is.null(Cfixed) & is.null(Dfixed)){
      
      # put the scale in Mode D
      adg <- colSums(Gnew^2)
      Gnew <- Gnew%*%(diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bnew^2)
      Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
      cdg <- colMeans(Cnew^2)
      Cnew <- Cnew%*%(diag(nfac)*(cdg^-0.5))
      Dnew <- Dnew%*%(diag(nfac)*((adg*bdg*cdg)^0.5))
      
      # order according to sum-of-squares
      fordr <- order(colSums(Dnew^2),decreasing=TRUE)
      Gnew <- as.matrix(Gnew[,fordr])
      Bnew <- as.matrix(Bnew[,fordr])
      Cnew <- as.matrix(Cnew[,fordr])
      Dnew <- as.matrix(Dnew[,fordr])
      
    }
    
    ### GCV criterion
    ntotal <- sum(sapply(data,function(x) dim(x)[1]))
    Adf <- nfac * (ntotal - xdims[4]*(nfac+1)/2)
    if(!is.null(Gstruc)) GtG <- crossprod(Gstruc)
    Gdf <- ifelse(is.null(Gstruc),
                  ifelse(const[1]==1L, 0, nfac*(nfac-1)/2),
                  sum(GtG[lower.tri(GtG)]>0L))
    Bdf <- ifelse(is.null(Bstruc),
                  nfac * ifelse(const[2]==1L, xdims[2]-(nfac+1)/2, xdims[2]-1L),
                  sum(Bstruc) - nfac)
    Cdf <- ifelse(is.null(Cstruc),
                  nfac * ifelse(const[3]==1L, xdims[3]-(nfac+1)/2, xdims[3]-1L),
                  sum(Cstruc) - nfac)
    Ddf <- ifelse(is.null(Dstruc),
                  nfac * ifelse(const[4]==1L, xdims[4]-(nfac-1)/2, xdims[4]),
                  sum(Dstruc))
    edf <- c(Adf+Gdf,Bdf,Cdf,Ddf)
    pxdim <- ntotal * prod(xdims[2:3])
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C","D")
    pfac <- list(A=list(H=Rknew,G=Gnew),B=Bnew,C=Cnew,D=Dnew,Rsq=Rsq,GCV=GCV,
                 edf=edf,iter=iter,cflag=cflag,const=const)
    return(pfac)
    
  }