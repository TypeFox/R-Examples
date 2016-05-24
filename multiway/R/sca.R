sca <- 
  function(X,nfac,nstart=10,maxit=500,
           type=c("sca-p","sca-pf2","sca-ind","sca-ecp"),
           rotation=c("none","varimax","promax"),
           ctol=10^-4,parallel=FALSE,cl=NULL){
    # Simultaneous Component Analysis
    # via alternating least squares (ALS) or closed-form solution
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 16, 2016
    
    # check 'X' input
    if(is.array(X)){
      xdim <- dim(X)
      lxdim <- length(xdim)
      if(lxdim!=3L){stop("Input 'X' must be 3-way array or list of 3-way arrays")}
      if(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
      mylist <- vector("list",xdim[3])
      for(kk in 1:xdim[3]){mylist[[kk]] <- X[,,kk]}
      X <- mylist
      rm(mylist)
    } else if(is.list(X)){
      d1x <- dim(X[[1]])
      lxdim <- length(d1x) + 1L
      if( any(sapply(X,function(x) any(any(is.na(x)),any(is.nan(x)),any(is.infinite(x))))) ){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
      if(lxdim==3L){
        xdim <- rep(NA,3)
        xdim[2] <- d1x[2]
        xdim[3] <- length(X)
        if(sum((sapply(X,ncol)-xdim[2])^2)>0L){stop("Input 'X' must be list of matrices with same number of columns.")}
      } else{stop("Input 'X' must be list of 2-way arrays.")}
    } else{stop("Input 'X' must be an array or list.")}
    xcx <- sumsq(X)
    
    # check 'nfac' and 'nstart' inputs
    nfac <- as.integer(nfac[1])
    if(nfac<1L){stop("Input 'nfac' must be positive integer")}
    nstart <- as.integer(nstart[1])
    if(nstart<1L){stop("Input 'nstart' must be positive integer")}
    
    # check 'maxit' and 'ctol' inputs
    maxit <- as.integer(maxit[1])
    if(maxit<1L){stop("Input 'maxit' must be positive integer")}
    ctol <- as.numeric(ctol[1])
    if(ctol<0L){stop("Input 'ctol' must be positive numeric")}
    
    # check 'type' and 'rotation'
    type <- type[1]
    if(is.na(match(type,c("sca-p","sca-pf2","sca-ind","sca-ecp")))){stop("Input 'type' must be one of four specified options")}
    rotation <- rotation[1]
    if(is.na(match(rotation,c("none","varimax","promax")))){stop("Input 'rotation' must be one of three specified options")}
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?sca")
    }
    
    # determine type of model fit
    if(type[1]=="sca-p"){
      matdata <- X[[1]]
      for(kk in 2:xdim[3]){matdata <- rbind(matdata,X[[kk]])}
      mysvd <- svd(matdata,nu=nfac,nv=nfac)
      Dmat <- mysvd$u%*%(diag(nfac)*mysvd$d[1:nfac]/sqrt(xdim[2]))
      Bmat <- mysvd$v*sqrt(xdim[2])
      ssr <- sum(rowSums(Dmat^2))
      if(rotation=="varimax"){
        Vrot <- varimax(Bmat)
        Bmat <- Bmat%*%Vrot$rotmat
        Dmat <- Dmat%*%Vrot$rotmat
      } else if(rotation=="promax"){
        Prot <- promax(Bmat)
        Bmat <- Bmat%*%Prot$rotmat
        Dmat <- Dmat%*%t(solve(Prot$rotmat))
      } 
      Dmats <- vector("list",xdim[3])
      dfc <- sse <- 0
      Cmat <- matrix(0,xdim[3],nfac)
      for(kk in 1:xdim[3]){
        newdim <- dim(X[[kk]])[1]
        dinds <- 1:newdim + dfc
        Dmats[[kk]] <- as.matrix(Dmat[dinds,])
        Cmat[kk,] <- sqrt(colSums(Dmats[[kk]]^2)/newdim)
        sse <- sse + sumsq(X[[kk]]-tcrossprod(Dmats[[kk]],Bmat))
        dfc <- dfc + newdim
      }
      rm(Dmat)
      Rsq <- 1 - sse/xcx
      iter <- 1
      cflag <- 0
      Phimat <- NULL
      ntotal <- nrow(matdata)
      Adf <- ntotal - xdim[3]*(nfac+1)/2
      Gdf <- xdim[3]*(nfac+1)/2
      Bdf <- xdim[2]-1L
      Cdf <- 0
      edf <- nfac * c(Adf+Gdf,Bdf,Cdf)
      pxdim <- ntotal * xdim[2]
      GCV <- (sse/pxdim) / (1 - sum(edf)/pxdim)^2
    } else if(type[1]=="sca-pf2"){
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac2_3way",data=X,xcx=xcx,
                              maxit=maxit,ctol=ctol)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac2_3way(data=X,nfac=nfac,xcx=xcx,maxit=maxit,ctol=ctol)
        }
      }
      widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
      scamod <- pfaclist[[widx]]
      rm(pfaclist)
      Bmat <- scamod$B
      Cmat <- matrix(0,xdim[3],nfac)
      Dmats <- vector("list",xdim[3])
      gsqrt <- sqrt(colSums(scamod$A$G^2))
      for(kk in 1:xdim[3]){
        Dmats[[kk]] <- scamod$A$H[[kk]]%*%scamod$A$G%*%(diag(nfac)*scamod$C[kk,])
        Cmat[kk,] <- scamod$C[kk,]*gsqrt/sqrt(dim(scamod$A$H[[kk]])[1])
      }
      Rsq <- scamod$Rsq
      iter <- scamod$iter
      cflag <- scamod$cflag
      Phimat <- crossprod(scamod$A$G%*%(diag(nfac)/gsqrt))
      GCV <- scamod$GCV
      edf <- scamod$edf
    } else if(type[1]=="sca-ind"){
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac2_3way",data=X,xcx=xcx,
                              const=c(1L,0L,0L),maxit=maxit,ctol=ctol)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac2_3way(data=X,nfac=nfac,xcx=xcx,const=c(1L,0L,0L),maxit=maxit,ctol=ctol)
        }
      }
      widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
      scamod <- pfaclist[[widx]]
      rm(pfaclist)
      dg <- sqrt(colSums(scamod$A$G^2))
      Bmat <- scamod$B
      Cmat <- matrix(0,xdim[3],nfac)
      Dmats <- vector("list",xdim[3])
      for(kk in 1:xdim[3]){
        Dmats[[kk]] <- scamod$A$H[[kk]]%*%scamod$A$G%*%(diag(nfac)*scamod$C[kk,])
        Cmat[kk,] <- scamod$C[kk,]*(dg/sqrt(dim(scamod$A$H[[kk]])[1]))
      }
      Rsq <- scamod$Rsq
      iter <- scamod$iter
      cflag <- scamod$cflag
      Phimat <- diag(nfac)
      GCV <- scamod$GCV
      edf <- scamod$edf
    } else if(type[1]=="sca-ecp"){
      nks <- rep(0,xdim[3])
      for(kk in 1:xdim[3]){nks[kk] <- sqrt(dim(X[[kk]])[1])}
      Cfixed <- matrix(nks,xdim[3],nfac)
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac2_3way",data=X,xcx=xcx,
                              const=c(1L,0L,0L),maxit=maxit,ctol=ctol,Cfixed=Cfixed)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac2_3way(data=X,nfac=nfac,xcx=xcx,const=c(1L,0L,0L),maxit=maxit,ctol=ctol,Cfixed=Cfixed)
        }
      }
      widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
      scamod <- pfaclist[[widx]]
      rm(pfaclist)
      fordr <- order(colSums(scamod$B^2), decreasing = TRUE)
      scamod$A$G <- as.matrix(scamod$A$G[,fordr])
      scamod$B <- as.matrix(scamod$B[,fordr])
      scamod$C <- as.matrix(scamod$C[,fordr])
      Bmat <- scamod$B
      rotmat <- diag(nfac)
      if(rotation=="varimax"){
        rotmat <- varimax(Bmat)$rotmat
        Bmat <- Bmat%*%rotmat
      } else if(rotation=="promax"){
        rotmat <- promax(Bmat)$rotmat
        Bmat <- Bmat%*%rotmat
        rotmat <- t(solve(rotmat))
      } 
      Dmats <- vector("list",xdim[3])
      for(kk in 1:xdim[3]){
        Dmats[[kk]] <- scamod$A$H[[kk]]%*%scamod$A$G%*%(diag(nfac)*scamod$C[kk,])%*%rotmat
      }
      Cmat <- scamod$C/nks
      Rsq <- scamod$Rsq
      iter <- scamod$iter
      cflag <- scamod$cflag
      Phimat <- crossprod(scamod$A$G)
      ntotal <- sum(sapply(X,nrow))
      Adf <- ntotal - xdim[3]*(nfac+1)/2
      Gdf <- 0
      Bdf <- xdim[2]
      Cdf <- 0
      edf <- nfac * c(Adf+Gdf,Bdf,Cdf)
      pxdim <- ntotal * xdim[2]
      ssenew <- (1-scamod$Rsq)*xcx
      GCV <- (ssenew/pxdim) / (1 - sum(edf)/pxdim)^2
    }
    names(edf) <- c("A","B","C")
    scafit <- list(D=Dmats,B=Bmat,C=Cmat,Phi=Phimat,
                   Rsq=Rsq,GCV=GCV,edf=edf,iter=iter,cflag=cflag,
                   type=type,rotation=rotation)
    class(scafit) <- "sca"
    return(scafit)
    
  }