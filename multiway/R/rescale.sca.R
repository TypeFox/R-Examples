rescale.sca <-
  function(x, mode="B", newscale=1, ...){
    # Rescales Weights of fit SCA nmodel
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: August 19, 2015    
    
    # check mode
    mode <- mode[1]
    if(mode=="B"){
      absorb <- "C"
    } else if(mode=="C"){
      absorb <- "B"
    } else {
      stop("Incorrect input for 'mode'. Must set to 'B' or 'C' for SCA")
    }
    
    # check newscale
    nfac <- ncol(x$B)
    if(length(newscale)!=nfac) newscale <- rep(newscale[1],nfac)
    
    # rescale factors
    if(mode=="B"){
      
      Bscale <- sqrt(colMeans(x$B^2))
      svec <- newscale/Bscale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      if(nfac==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      x$C <- x$C %*% Smat
      x$D <- lapply(x$D, function(x) x %*% Smat)
      if(x$type=="sca-ecp") { x$Phi <- Smat %*% x$Phi %*% Smat }
      return(x)
      
    } else {
      
      Cscale <- sqrt(colMeans(x$C^2))
      svec <- newscale/Cscale
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      x$D <- lapply(x$D, function(x) x %*% Smat)
      if(x$type=="sca-ecp") { x$Phi <- Smat %*% x$Phi %*% Smat }
      if(nfac==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      x$B <- x$B %*% Smat
      return(x)
      
    }
    
  }