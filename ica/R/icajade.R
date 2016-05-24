icajade <-
  function(X,nc,center=TRUE,maxit=100,tol=1e-6,Rmat=diag(nc)){
    ###### Joint Approximate Diagonalization of Eigenmatrices (JADE)
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: August 23, 2015
    
    ### initial checks
    X <- as.matrix(X)
    nobs <- nrow(X)
    nvar <- ncol(X)
    nc <- as.integer(nc[1])
    if(nc<1){ stop("Must set nc>=1 component.") }
    maxit <- as.integer(maxit[1])
    if(maxit<1){ stop("Must set maxit>=1 iteration.") }
    tol <- tol[1]
    if(tol<=0){ stop("Must set ctol>0.") }
    if(nc>min(nobs,nvar)){ stop("Too many components. Set nc<=min(dim(X)).") }
    if(nrow(Rmat)!=nc | ncol(Rmat)!=nc){ stop("Input 'Rmat' must be nc-by-nc rotation matrix.") }
    
    ### center and whiten
    if(center) X <- scale(X,scale=FALSE)
    xeig <- eigen(crossprod(X)/nobs,symmetric=TRUE)
    nze <- sum(xeig$val>xeig$val[1]*.Machine$double.eps)
    if(nze<nc){
      warning("Numerical rank of X is less than requested number of components (nc).\n  Number of components has been redefined as the numerical rank of X.")
      nc <- nze
      Rmat <- diag(nc)
    }
    Dmat <- sdiag(sqrt(xeig$val[1:nc]))
    Mprt <- tcrossprod(Dmat,xeig$vec[,1:nc])
    diag(Dmat) <- 1/diag(Dmat)
    Pmat <- xeig$vec[,1:nc]%*%Dmat
    Xw <- X%*%Pmat   # whitened data
    
    ### check if nc=1
    if(nc==1L){
      return(list(S=Xw,M=Mprt,W=t(Pmat),Y=Xw,Q=t(Pmat),R=matrix(1),
                  vafs=(sum(Mprt^2)*nobs)/sum(X^2),iter=NA))
    }
    
    ### basis eigenmatrices (using Jean-Francois Cardoso's symmetry trick)
    ncstar <- nc*(nc+1)/2
    idmat <- diag(nc)
    emats <- matrix(0,nc,nc*ncstar)
    crng <- 1:nc
    for(i in 1:nc){
      Xi <- Xw[,i]
      Qij <- crossprod(matrix((Xi^2)/nobs,nobs,nc)*Xw,Xw) - idmat - 2*tcrossprod(idmat[,i],idmat[,i])
      emats[,crng] <- Qij
      crng <- crng + nc
      if(i>1){
        for(j in 1:(i-1)){
          Xj <- Xw[,j]
          Qij <- crossprod(matrix(Xi*Xj/nobs,nobs,nc)*Xw,Xw) - tcrossprod(idmat[,i],idmat[,j]) - tcrossprod(idmat[,j],idmat[,i])
          emats[,crng] <- sqrt(2)*Qij
          crng <- crng + nc
        } # end if(i>1)
      } # end for(j in 1:(i-1))
    } # end for(i in 1:nc)
    
    ### iterative rotation
    npairs <- nc*(nc-1)/2
    thetas <- rep(1,npairs)
    iter <- 0
    vtol <- 1
    while(vtol>tol && iter<maxit){
      # sweep through angle pairs
      for(p in 1:(nc-1)){
        for(q in (p+1):nc){
          # Givens angle
          ip <- seq(p,nc*ncstar,by=nc)
          iq <- seq(q,nc*ncstar,by=nc)
          gp <- rbind(emats[p,ip]-emats[q,iq],emats[p,iq]+emats[q,ip])
          gg <- tcrossprod(gp)
          ton <- gg[1,1]-gg[2,2]
          toff <- gg[1,2]+gg[2,1]
          theta <- 0.5*atan2(toff,ton+sqrt(ton^2+toff^2))
          thetas[nc*(p-1)-p*(p-1)/2+q-p] <- theta
          # Givens rotation
          cc <- cos(theta)
          ss <- sin(theta)
          gmat <- rbind(c(cc,-ss),c(ss,cc))
          pair <- c(p,q)
          Rmat[,pair] <- Rmat[,pair]%*%gmat
          emats[pair,] <- crossprod(gmat,emats[pair,])
          emats[,c(ip,iq)] <- cbind(cc*emats[,ip]+ss*emats[,iq],-ss*emats[,ip]+cc*emats[,iq])
        }
      }
      # check for convergence
      vtol <- max(abs(thetas))
      iter <- iter + 1
    } # end while(vtol>tol && iter<maxit)
    
    ### sort according to vafs
    M <- crossprod(Rmat,Mprt)
    vafs <- rowSums(M^2)
    ix <- sort(vafs,decreasing=TRUE,index.return=TRUE)$ix
    M <- M[ix,]
    Rmat <- Rmat[,ix]
    vafs <- (vafs[ix]*nobs)/sum(X^2)
    
    return(list(S=Xw%*%Rmat,M=t(M),W=t(Pmat%*%Rmat),Y=Xw,
                Q=t(Pmat),R=Rmat,vafs=vafs,iter=iter))
    
  }