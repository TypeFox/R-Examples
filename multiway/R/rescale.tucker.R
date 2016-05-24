rescale.tucker <-
  function(x, mode="A", newscale=1, ...){
    # Rescales Weights of fit Tucker model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: August 19, 2015    
    
    # check mode and dimensions
    mode <- mode[1]
    mydim <- c(ncol(x$A),ncol(x$B),ncol(x$C))
    apdim <- 1:3
    if(is.null(x$D)){
      if(!any(mode==c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C' for 3-way Tucker")
    } else {
      mydim <- c(mydim,ncol(x$D))
      apdim <- c(apdim,4)
      if(!any(mode==c("A","B","C","D"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', 'C', or 'D' for 4-way Parafac.")
    }
    
    # rescale factors
    if(mode=="A"){
      
      if(length(newscale)!=mydim[1]) newscale <- rep(newscale[1],mydim[1])
      Ascale <- sqrt(colMeans(x$A^2))
      svec <- newscale/Ascale
      if(mydim[1]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$A <- x$A %*% Smat
      if(mydim[1]==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      Gmat <- matrix(x$G, mydim[1], prod(mydim[-1]))
      x$G <- array(Smat %*% Gmat, dim=mydim)
      return(x)
      
    } else if(mode=="B"){
      
      if(length(newscale)!=mydim[2]) newscale <- rep(newscale[1],mydim[2])
      permvec <- c(apdim[2],apdim[-2])
      Bscale <- sqrt(colMeans(x$B^2))
      svec <- newscale/Bscale
      if(mydim[2]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$B <- x$B %*% Smat
      if(mydim[2]==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      Gmat <- matrix(aperm(x$G, permvec), mydim[2], prod(mydim[-2]))
      x$G <- aperm(array(Smat %*% Gmat, dim=c(mydim[2],mydim[-2])), sort(permvec,index=T)$ix)
      return(x)
      
    } else if(mode=="C"){
      
      if(length(newscale)!=mydim[3]) newscale <- rep(newscale[1],mydim[3])
      permvec <- c(apdim[3],apdim[-3])
      Cscale <- sqrt(colMeans(x$C^2))
      svec <- newscale/Cscale
      if(mydim[3]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$C <- x$C %*% Smat
      if(mydim[3]==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      Gmat <- matrix(aperm(x$G, permvec), mydim[3], prod(mydim[-3]))
      x$G <- aperm(array(Smat %*% Gmat, dim=c(mydim[3],mydim[-3])), sort(permvec,index=T)$ix)
      return(x)
      
    } else if(mode=="D"){
      
      if(length(newscale)!=mydim[4]) newscale <- rep(newscale[1],mydim[4])
      permvec <- c(apdim[4],apdim[-4])
      Dscale <- sqrt(colMeans(x$D^2))
      svec <- newscale/Dscale
      if(mydim[4]==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
      x$D <- x$D %*% Smat
      if(mydim[4]==1L) { Smat <- matrix(1/svec) } else { Smat <- diag(1/svec) }
      Gmat <- matrix(aperm(x$G, permvec), mydim[4], prod(mydim[-4]))
      x$G <- aperm(array(Smat %*% Gmat, dim=c(mydim[4],mydim[-4])), sort(permvec,index=T)$ix)
      return(x)
      
    }
    
  }