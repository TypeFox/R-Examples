smartssp <-
  function(Qmats,lambdas,Kty,Jty,KtK,KtJ,
           JtJ,nknots,ndpts,alpha,yty,nbf){
    # theta smart start for flexible ssa functions
    
    # set thetas so subspaces have equal influence
    nfs <- dim(Qmats)[2]/nknots
    sdq <- rep(NA,nfs)
    for(jj in 1:nfs){
      jind <- ((jj-1)*nknots+1):(jj*nknots)
      sdq[jj] <- sum(diag(Qmats[,jind]))
    }
    thetas <- 1/sdq
    
    # fit one smoothing parameter model with given thetas
    thvec <- thetas
    chat <- (lamcoef(lambdas,thvec,Kty,Jty,KtK,KtJ,JtJ,
                     Qmats,nknots,ndpts,alpha,yty,nbf))[[1]][(nbf+1):(nknots+nbf)]
    
    # solve for optimal thetas
    for(jj in 1:nfs){
      jind <- ((jj-1)*nknots+1):(jj*nknots)
      thetas[jj] <- (thetas[jj]^2)*crossprod(pdsXty(Qmats[,jind],chat))
    }
    
    thetas
    
  }
