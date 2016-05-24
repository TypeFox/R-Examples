resign.sca <-
  function(x, mode="B", newsign=1, ...){
    # Resigns Weights of fit SCA nmodel
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
    
    # check newsign
    nfac <- ncol(x$B)
    newsign <- sign(newsign)
    if(length(newsign)!=nfac) newsign <- rep(newsign[1],nfac)
    
    # resign factors
    if(mode=="B"){
      
      Bsign <- sign(colMeans(x$B^3))
      svec <- newsign*Bsign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      x$C <- x$C %*% Smat
      x$D <- lapply(x$D, function(x) x %*% Smat)
      if(x$type=="sca-ecp") { x$Phi <- Smat %*% x$Phi %*% Smat }
      return(x)
      
    } else {
      
      Csign <- sign(colMeans(x$C^3))
      svec <- newsign*Csign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      x$D <- lapply(x$D, function(x) x %*% Smat)
      if(x$type=="sca-ecp") { x$Phi <- Smat %*% x$Phi %*% Smat }
      x$B <- x$B %*% Smat
      return(x)
      
    }
    
  }