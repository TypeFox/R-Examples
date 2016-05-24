tucker <- 
  function(X,nfac,nstart=10,Afixed=NULL,
           Bfixed=NULL,Cfixed=NULL,Dfixed=NULL,
           Bstart=NULL,Cstart=NULL,Dstart=NULL,
           maxit=500,ctol=10^-4,parallel=FALSE,
           cl=NULL,output=c("best","all")){
    # 3-way or 4-way Tucker model
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
    nfac <- as.integer(nfac)
    if(length(nfac)!=lxdim){stop(paste("Input 'nfac' must be ",lxdim," element vector specifying number of factors for each mode"))}
    if(any(nfac<1L)){stop("Input 'nfac' must contain positive integers")}
    nstart <- as.integer(nstart[1])
    if(nstart<1L){stop("Input 'nstart' must be positive integer")}
    
    # check 'maxit' and 'ctol' inputs
    maxit <- as.integer(maxit[1])
    if(maxit<1L){stop("Input 'maxit' must be positive integer")}
    ctol <- as.numeric(ctol[1])
    if(ctol<0L){stop("Input 'ctol' must be positive numeric")}
    
    # check 'Afixed', 'Bfixed', and 'Cfixed' inputs
    if(!is.null(Afixed)){
      Afixed <- as.matrix(Afixed)
      if(nrow(Afixed)!=xdim[1]){stop("Input 'Afixed' must have the same number of rows as dim(X)[1]")}
      if(ncol(Afixed)!=nfac){stop("Input 'Afixed' must have 'nfac' columns")}
    }
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
    
    # check 'parallel' and 'cl' inputs
    if(parallel && !any(class(cl)=="cluster")){
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?tucker")
    }
    
    # tucker fitting
    if(lxdim==3L){
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist <- lapply(nstartlist,function(x) {x <- nfac})
        tucklist <- parLapply(cl=cl,X=nstartlist,fun="tucker_3way",data=X,xcx=xcx,
                              maxit=maxit,ctol=ctol,Afixed=Afixed,Bfixed=Bfixed,Cfixed=Cfixed,
                              Bstart=Bstart,Cstart=Cstart)
      } else {
        tucklist <- vector("list",nstart)
        for(j in 1:nstart){
          tucklist[[j]] <- tucker_3way(data=X,nfac=nfac,xcx=xcx,maxit=maxit,ctol=ctol,
                                       Afixed=Afixed,Bfixed=Bfixed,Cfixed=Cfixed,
                                       Bstart=Bstart,Cstart=Cstart)
        }
      } # end if(parallel)
    } else if(lxdim==4L){
      # check 'Dfixed' and 'Dstart' inputs
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
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist <- lapply(nstartlist,function(x) {x <- nfac})
        tucklist <- parLapply(cl=cl,X=nstartlist,fun="tucker_4way",data=X,xcx=xcx,
                              maxit=maxit,ctol=ctol,Afixed=Afixed,Bfixed=Bfixed,Cfixed=Cfixed,
                              Dfixed=Dfixed,Bstart=Bstart,Cstart=Cstart,Dstart=Dstart)
      } else {
        tucklist <- vector("list",nstart)
        for(j in 1:nstart){
          tucklist[[j]] <- tucker_4way(data=X,nfac=nfac,xcx=xcx,maxit=maxit,ctol=ctol,
                                       Afixed=Afixed,Bfixed=Bfixed,Cfixed=Cfixed,Dfixed=Dfixed,
                                       Bstart=Bstart,Cstart=Cstart,Dstart=Dstart)
        }
      } # end if(parallel)
    } # end if(lxdim==3L) 
    
    # output results
    if(output[1]=="best"){
      widx <- which.max(sapply(tucklist,function(x) x$Rsq))
      tuck <- tucklist[[widx]]
      class(tuck) <- "tucker"
      return(tuck)
    } else {
      tucklist <- lapply(tucklist, function(tuck) {
        class(tuck) <- "tucker"
        tuck
        })
      return(tucklist)
    }
    
  } # end tucker.R