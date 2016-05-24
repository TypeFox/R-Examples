parafac_4way <-
  function(data,nfac,xcx=sumsq(data),const=rep(0L,4),
           maxit=500,ctol=10^-4,Bfixed=NULL,Cfixed=NULL,
           Dfixed=NULL,Bstart=NULL,Cstart=NULL,Dstart=NULL,
           Bstruc=NULL,Cstruc=NULL,Dstruc=NULL){
    # 4-way Parallel Factor Analysis (Parafac)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: January 27, 2016
    
    ### initialize Khatri-Rao product matrices
    xdims <- dim(data)
    if(is.null(Dfixed)){CBkrA <- matrix(0,xdims[1]*xdims[2]*xdims[3],nfac)}
    if(is.null(Cfixed)){DBkrA <- matrix(0,xdims[1]*xdims[2]*xdims[4],nfac)}
    if(is.null(Bfixed)){DCkrA <- matrix(0,xdims[1]*xdims[3]*xdims[4],nfac)}
    DCkrB <- matrix(0,xdims[2]*xdims[3]*xdims[4],nfac)
    
    ### initialize reshaped data matrices
    Xa <- matrix(data,xdims[1],xdims[2]*xdims[3]*xdims[4])
    if(is.null(Bfixed)){Xb <- matrix(aperm(data,perm=c(2,1,3,4)),xdims[2],xdims[1]*xdims[3]*xdims[4])}
    if(is.null(Cfixed)){Xc <- matrix(aperm(data,perm=c(3,1,2,4)),xdims[3],xdims[1]*xdims[2]*xdims[4])}
    if(is.null(Dfixed)){Xd <- matrix(aperm(data,perm=c(4,1,2,3)),xdims[4],xdims[1]*xdims[2]*xdims[3])}
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
    if(is.null(Dfixed)){
      if(!is.null(Dstart)){
        Dold <- Dstart
      } else if(const[4]==0L){
        Dold <- matrix(rnorm(xdims[4]*nfac),xdims[4],nfac)
      } else if(const[4]==1L){
        Dold <- svd(matrix(rnorm(xdims[4]*nfac),xdims[4],nfac),nu=nfac,nv=0)$u
      } else if(const[4]==2L){
        Dold <- Dnew <- matrix(runif(xdims[4]*nfac),xdims[4],nfac)
      }
      if(!is.null(Dstruc)) Dold <- Dnew <- Dold * Dstruc
    } else {Dold=Dnew=Dfixed}
    
    ### iterative update of matrices
    vtol <- sseold <- xcx
    iter <- 0
    cflag <- NA
    while(vtol>ctol && iter<maxit) {
      
      # Step 1: update mode A weights
      for(u in 1:nfac){DCkrB[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Bold[,u]))}
      if(const[1]==0L){
        #Anew <- Xa%*%DCkrB%*%smpower(crossprod(DCkrB),-1)
        Anew <- Xa%*%DCkrB%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Bold),-1)
      } else if(const[1]==1L) {
        Zmat <- Xa%*%DCkrB
        Anew <- Zmat%*%smpower(crossprod(Zmat),-0.5)
      } else if(const[1]==2L) {
        cpmat <- crossprod(DCkrB)
        for(ii in 1:xdims[1]){Anew[ii,] <- fnnls(cpmat,crossprod(DCkrB,Xa[ii,]))}
        if(any(colSums(Anew)==0)){
          Anew <- Aold
          vtol <- 0
          cflag <- 2
        }
      }
      
      # Step 2: update mode B weights
      if(is.null(Bfixed)){
        for(u in 1:nfac){DCkrA[,u] <- kronecker(Dold[,u],kronecker(Cold[,u],Anew[,u]))}
        if(is.null(Bstruc)){
          if(const[2]==0L){
            #Bnew <- Xb%*%DCkrA%*%smpower(crossprod(DCkrA),-1)
            Bnew <- Xb%*%DCkrA%*%smpower(crossprod(Dold)*crossprod(Cold)*crossprod(Anew),-1)
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
      
      # Step 3: update mode C weights
      if(is.null(Cfixed)){
        for(u in 1:nfac){DBkrA[,u] <- kronecker(Dold[,u],kronecker(Bnew[,u],Anew[,u]))}
        if(is.null(Cstruc)){
          if(const[3]==0L){
            #Cnew <- Xc%*%DBkrA%*%smpower(crossprod(DBkrA),-1)
            Cnew <- Xc%*%DBkrA%*%smpower(crossprod(Dold)*crossprod(Bnew)*crossprod(Anew),-1)
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
      
      # Step 4: update mode D weights
      if(is.null(Dfixed)){
        for(u in 1:nfac){CBkrA[,u] <- kronecker(Cnew[,u],kronecker(Bnew[,u],Anew[,u]))}
        if(is.null(Dstruc)){
          if(const[4]==0L){
            #Dnew <- Xd%*%CBkrA%*%smpower(crossprod(CBkrA),-1)
            Dnew <- Xd%*%CBkrA%*%smpower(crossprod(Cnew)*crossprod(Bnew)*crossprod(Anew),-1)
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
      
      # Step 5: check for convergence
      for(u in 1:nfac){DCkrB[,u] <- kronecker(Dnew[,u],kronecker(Cnew[,u],Bnew[,u]))}
      ssenew <- sum((Xa-tcrossprod(Anew,DCkrB))^2)
      #vtol <- (sseold-ssenew)/sseold
      vtol <- (sseold - ssenew) / xcx
      Aold <- Anew
      Bold <- Bnew
      Cold <- Cnew
      Dold <- Dnew
      sseold <- ssenew
      iter <- iter + 1
      
    }  # end while(vtol>ctol && iter<maxit)
    
    ### scale and order solution
    if(is.null(Bfixed) & is.null(Cfixed) & is.null(Dfixed)){
      
      # put the scale in Mode D
      adg <- colMeans(Anew^2)
      Anew <- Anew%*%(diag(nfac)*(adg^-0.5))
      bdg <- colMeans(Bnew^2)
      Bnew <- Bnew%*%(diag(nfac)*(bdg^-0.5))
      cdg <- colMeans(Cnew^2)
      Cnew <- Cnew%*%(diag(nfac)*(cdg^-0.5))
      Dnew <- Dnew%*%(diag(nfac)*((adg*bdg*cdg)^0.5))
      
      # order according to sum-of-squares
      fordr <- order(colSums(Dnew^2),decreasing=TRUE)
      Anew <- as.matrix(Anew[,fordr])
      Bnew <- as.matrix(Bnew[,fordr])
      Cnew <- as.matrix(Cnew[,fordr])
      Dnew <- as.matrix(Dnew[,fordr])
      
    }
    
    ### GCV criterion
    Adf <- nfac * ifelse(const[1]==1L, xdims[1]-(nfac+1)/2, xdims[1]-1L)
    Bdf <- ifelse(is.null(Bstruc),
                  nfac * ifelse(const[2]==1L, xdims[2]-(nfac+1)/2, xdims[2]-1L),
                  sum(Bstruc) - nfac)
    Cdf <- ifelse(is.null(Cstruc),
                  nfac * ifelse(const[3]==1L, xdims[3]-(nfac+1)/2, xdims[3]-1L),
                  sum(Cstruc) - nfac)
    Ddf <- ifelse(is.null(Dstruc),
                  nfac * ifelse(const[4]==1L, xdims[4]-(nfac-1)/2, xdims[4]),
                  sum(Dstruc))
    edf <- c(Adf,Bdf,Cdf,Ddf)
    pxdim <- prod(xdims)
    GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    
    ### collect results
    Rsq <- 1 - ssenew/xcx
    if(is.na(cflag)){
      if(vtol<=ctol){cflag <- 0} else {cflag <- 1}
    }
    names(edf) <- c("A","B","C","D")
    pfac <- list(A=Anew,B=Bnew,C=Cnew,D=Dnew,Rsq=Rsq,GCV=GCV,edf=edf,
                 iter=iter,cflag=cflag,const=const)
    return(pfac)
    
  }