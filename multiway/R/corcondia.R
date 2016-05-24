corcondia <- 
  function(X,object,divisor=c("nfac","core")){
    # Core Consistency Diagnostic
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: October 12, 2015
    
    clob <- class(object)
    if(clob=="parafac"){
      
      nfac <- ncol(object$A)
      Ai <- mpinv(object$A)
      Bi <- mpinv(object$B)
      Ci <- mpinv(object$C)
      imat <- kronecker(Ci,kronecker(Bi,Ai))
      if(is.null(object$D)){
        g <- imat %*% c(X)
        G <- array(g, dim=rep(nfac,3))
        super <- array(0, dim=rep(nfac,3))
        for(k in 1:nfac) super[k,k,k] <- 1
      } else {
        Di <- mpinv(object$D)
        imat <- kronecker(Di,imat)
        g <- imat %*% c(X)
        G <- array(g, dim=rep(nfac,4))
        super <- array(0, dim=rep(nfac,4))
        for(k in 1:nfac) super[k,k,k,k] <- 1
      }
      
    } else if(clob=="parafac2"){
      
      if(is.array(X)){
        xdim <- dim(X)
        lxdim <- length(xdim)
        if(lxdim==3L){
          Xlist <- vector("list",nrow(object$C))
          for(k in 1:nrow(object$C)) Xlist[[k]] <- X[,,k]
        } else {
          Xlist <- vector("list",nrow(object$D))
          for(k in 1:nrow(object$D)) Xlist[[k]] <- X[,,,k]
        }
        X <- Xlist
        rm(Xlist)
      }
      
      nfac <- ncol(object$A$G)
      Ai <- mpinv(object$A$G)
      Bi <- mpinv(object$B)
      Ci <- mpinv(object$C)
      imat <- kronecker(Ci,kronecker(Bi,Ai))
      
      if(is.null(object$D)){
        for(k in 1:nrow(object$C)) X[[k]] <- crossprod(object$A$H[[k]], X[[k]])
        g <- imat %*% unlist(X)
        G <- array(g, dim=rep(nfac,3))
        super <- array(0, dim=rep(nfac,3))
        for(k in 1:nfac) super[k,k,k] <- 1
      } else {
        Di <- mpinv(object$D)
        imat <- kronecker(Di,imat)
        for(k in 1:nrow(object$D)) X[[k]] <- crossprod(object$A$H[[k]], matrix(X[[k]],nrow=dim(X[[k]])[1]))
        g <- imat %*% unlist(X)
        G <- array(g, dim=rep(nfac,4))
        super <- array(0, dim=rep(nfac,4))
        for(k in 1:nfac) super[k,k,k,k] <- 1
      }
      
    } else {
      stop("Input 'object' must be object of class 'parafac' or 'parafac2'")
    }
    
    if(divisor[1]=="nfac"){
      corcon <- 100 * (1 - sum((G-super)^2)/nfac)
    } else {
      corcon <- 100 * (1 - sum((G-super)^2)/sum(G^2))
    }
    
    return(corcon)
  
}