parafac <- 
  function(X,nfac,nstart=10,const=NULL,
           Bfixed=NULL,Cfixed=NULL,Dfixed=NULL,
           Bstart=NULL,Cstart=NULL,Dstart=NULL,
           Bstruc=NULL,Cstruc=NULL,Dstruc=NULL,
           maxit=500,ctol=10^-4,parallel=FALSE,
           cl=NULL,output=c("best","all")){
    # 3-way or 4-way Parallel Factor Analysis (Parafac)
    # via alternating least squares (ALS) with optional constraints
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 4, 2016
    
    # check 'X' input
    xdim <- dim(X)
    lxdim <- length(xdim)
    if(lxdim<3L | lxdim>4L){stop("Input 'X' must be 3-way or 4-way array")}
    if(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
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
    
    # check 'const' input
    if(is.null(const)){
      const <- rep(0L,lxdim)
    } else {
      const <- as.integer(const)
      if(length(const)!=lxdim){stop(paste("Input 'const' must be ",lxdim," element vector specifying constraint for each mode"))}
      if(min(const)<0L | max(const)>2L){stop("'const[j]' must be 0 (unconstrained), 1 (orthogonal), or 2 (non-negative)")}
    } # end if(is.null(const))
    
    # check 'Bfixed' and 'Cfixed' inputs
    if(!is.null(Bfixed)){
      Bfixed <- as.matrix(Bfixed)
      if(nrow(Bfixed)!=xdim[2]){stop("Input 'Bfixed' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bfixed)!=nfac){stop("Input 'Bfixed' must have 'nfac' columns")}
    }
    if(!is.null(Cfixed)){
      Cfixed <- as.matrix(Cfixed)
      if(nrow(Cfixed)!=xdim[3]){stop("Input 'Cfixed' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cfixed)!=nfac){stop("Input 'Cfixed' must have 'nfac' columns")}
    }
    
    # check 'Bstart' and 'Cstart' inputs
    if(!is.null(Bstart)){
      Bstart <- as.matrix(Bstart)
      if(nrow(Bstart)!=xdim[2]){stop("Input 'Bstart' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bstart)!=nfac){stop("Input 'Bstart' must have 'nfac' columns")}
    }
    if(!is.null(Cstart)){
      Cstart <- as.matrix(Cstart)
      if(nrow(Cstart)!=xdim[3]){stop("Input 'Cstart' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cstart)!=nfac){stop("Input 'Cstart' must have 'nfac' columns")}
    }
    
    # check 'Bstruc' and 'Cstruc' inputs
    if(!is.null(Bstruc)){
      Bstruc <- as.matrix(Bstruc)
      if(nrow(Bstruc)!=xdim[2]){stop("Input 'Bstruc' must have the same number of rows as dim(X)[2]")}
      if(ncol(Bstruc)!=nfac){stop("Input 'Bstruc' must have 'nfac' columns")}
      if(any("logical"!=c(apply(Bstruc,1:2,class)))){stop("Input 'Bstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
    }
    if(!is.null(Cstruc)){
      Cstruc <- as.matrix(Cstruc)
      if(nrow(Cstruc)!=xdim[3]){stop("Input 'Cstruc' must have the same number of rows as dim(X)[3]")}
      if(ncol(Cstruc)!=nfac){stop("Input 'Cstruc' must have 'nfac' columns")}
      if(any("logical"!=c(apply(Cstruc,1:2,class)))){stop("Input 'Cstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
    }
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?parafac")
    }
    
    # parafac fitting
    if(lxdim==3L){
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac_3way",data=X,xcx=xcx,
                              const=const,maxit=maxit,ctol=ctol,Bfixed=Bfixed,Cfixed=Cfixed,
                              Bstart=Bstart,Cstart=Cstart,Bstruc=Bstruc,Cstruc=Cstruc)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac_3way(data=X,nfac=nfac,xcx=xcx,const=const,maxit=maxit,ctol=ctol,
                                        Bfixed=Bfixed,Cfixed=Cfixed,Bstart=Bstart,Cstart=Cstart,
                                        Bstruc=Bstruc,Cstruc=Cstruc)
        }
      } # end if(parallel)
    } else if(lxdim==4L){
      # check 'Dfixed', 'Dstart', and 'Dstruc' inputs
      if(!is.null(Dfixed)){
        Dfixed <- as.matrix(Dfixed)
        if(nrow(Dfixed)!=xdim[4]){stop("Input 'Dfixed' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dfixed)!=nfac){stop("Input 'Dfixed' must have 'nfac' columns")}
      }
      if(!is.null(Dstart)){
        Dstart <- as.matrix(Dstart)
        if(nrow(Dstart)!=xdim[4]){stop("Input 'Dstart' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dstart)!=nfac){stop("Input 'Dstart' must have 'nfac' columns")}
      }
      if(!is.null(Dstruc)){
        Dstruc <- as.matrix(Dstruc)
        if(nrow(Dstruc)!=xdim[4]){stop("Input 'Dstruc' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dstruc)!=nfac){stop("Input 'Dstruc' must have 'nfac' columns")}
        if(any("logical"!=c(apply(Dstruc,1:2,class)))){stop("Input 'Dstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
      }
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac_4way",data=X,xcx=xcx,
                              const=const,maxit=maxit,ctol=ctol,Bfixed=Bfixed,Cfixed=Cfixed,
                              Dfixed=Dfixed,Bstart=Bstart,Cstart=Cstart,Dstart=Dstart,
                              Bstruc=Bstruc,Cstruc=Cstruc,Dstruc=Dstruc)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac_4way(data=X,nfac=nfac,xcx=xcx,const=const,maxit=maxit,
                                        ctol=ctol,Bfixed=Bfixed,Cfixed=Cfixed,Dfixed=Dfixed,
                                        Bstart=Bstart,Cstart=Cstart,Dstart=Dstart,
                                        Bstruc=Bstruc,Cstruc=Cstruc,Dstruc=Dstruc)
        }
      } # end if(parallel)
    } # end if(lxdim==3L)
    
    # output results
    if(output[1]=="best"){
      widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
      pfac <- pfaclist[[widx]]
      class(pfac) <- "parafac"
      return(pfac)
    } else {
      pfaclist <- lapply(pfaclist, function(pfac) {
        class(pfac) <- "parafac"
        pfac
        })
      return(pfaclist)
    }
    
  } # end parafac.R