resign.parafac <-
  function(x, mode="A", newsign=1, absorb="C", ...){
    # Resigns Weights of fit Parafac model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: August 19, 2015    
    
    # check mode and absorb
    mode <- mode[1]
    absorb <- absorb[1]
    if(mode==absorb) stop("Inputs 'mode' and 'absorb' must be different.")
    if(is.null(x$D)){
      if(!any(mode==c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C' for 3-way Parafac.")
      if(!any(absorb==c("A","B","C"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', or 'C' for 3-way Parafac.")
    } else {
      if(!any(mode==c("A","B","C","D"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac.")
      if(!any(absorb==c("A","B","C","D"))) stop("Incorrect input for 'absorb'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac.")
    }
    
    # check newsign
    nfac <- ncol(x$A)
    newsign <- sign(newsign)
    if(length(newsign)!=nfac) newsign <- rep(newsign[1],nfac)
    
    # resign factors
    if(mode=="A"){
      
      Asign <- sign(colMeans(x$A^3))
      svec <- newsign*Asign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$A <- x$A %*% Smat
      if(absorb=="B") {
        x$B <- x$B %*% Smat
      } else if(absorb=="C"){
        x$C <- x$C %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="B"){
      
      Bsign <- sign(colMeans(x$B^3))
      svec <- newsign*Bsign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      if(absorb=="A") {
        x$A <- x$A %*% Smat
      } else if(absorb=="C"){
        x$C <- x$C %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="C"){
      
      Csign <- sign(colMeans(x$C^3))
      svec <- newsign*Csign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      if(absorb=="A") {
        x$A <- x$A %*% Smat
      } else if(absorb=="B"){
        x$B <- x$B %*% Smat
      } else {
        x$D <- x$D %*% Smat
      }
      return(x)
      
    } else if(mode=="D"){
      
      Dsign <- sign(colMeans(x$D^3))
      svec <- newsign*Dsign
      if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$D <- x$D %*% Smat
      if(absorb=="A") {
        x$A <- x$A %*% Smat
      } else if(absorb=="B"){
        x$B <- x$B %*% Smat
      } else {
        x$C <- x$C %*% Smat
      }
      return(x)
      
    }
    
  }