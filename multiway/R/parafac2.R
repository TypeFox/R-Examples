parafac2 <- 
  function(X,nfac,nstart=10,const=NULL,
           Gfixed=NULL,Bfixed=NULL,Cfixed=NULL,Dfixed=NULL,
           Gstart=NULL,Bstart=NULL,Cstart=NULL,Dstart=NULL,
           Gstruc=NULL,Bstruc=NULL,Cstruc=NULL,Dstruc=NULL,
           maxit=500,ctol=10^-4,parallel=FALSE,cl=NULL,output=c("best","all")){
    # 3-way or 4-way Parallel Factor Analysis2 (Parafac2)
    # via alternating least squares (ALS) with optional constraints
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: February 4, 2016
    
    # check 'X' input
    if(is.array(X)){
      xdim <- dim(X)
      lxdim <- length(xdim)
      if(lxdim<3L | lxdim>4L){stop("Input 'X' must be 3-way or 4-way array")}
      if(any(is.na(X)) | any(is.nan(X)) | any(is.infinite(X))){stop("Input 'X' cannot contain NA, NaN, or Inf values")}
      if(lxdim==3L){
        mylist <- vector("list",xdim[3])
        for(kk in 1:xdim[3]){mylist[[kk]] <- X[,,kk]}
      } else {
        mylist <- vector("list",xdim[4])
        for(kk in 1:xdim[4]){mylist[[kk]] <- X[,,,kk]}
      }
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
      } else if(lxdim==4L){
        xdim <- rep(NA,4)
        xdim[2] <- d1x[2]
        xdim[3] <- d1x[3]
        xdim[4] <- length(X)
        if(sum((sapply(X,dim)[2,]-xdim[2])^2)>0L){stop("Input 'X' must be list of arrays with same number of columns.")}
        if(sum((sapply(X,dim)[3,]-xdim[3])^2)>0L){stop("Input 'X' must be list of arrays with same number of slabs.")}
      } else{stop("Input 'X' must be list of 2-way or 3-way arrays.")}
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
    
    # check 'const' input
    if(is.null(const)){
      const <- rep(0L,lxdim)
    } else {
      const <- as.integer(const)
      if(length(const)!=lxdim){stop(paste("Input 'const' must be ",lxdim," element vector specifying constraint for each mode"))}
      if(min(const)<0L | max(const)>2L){stop("'const[j]' must be 0 (unconstrained), 1 (orthogonal), or 2 (non-negative)")}
      if(const[1]==2L){stop("you cannot apply nonnegativity on Mode A")}
    }
    
    # check 'Gfixed' and 'Bfixed' and 'Cfixed' inputs
    if(!is.null(Gfixed)){
      Gfixed <- as.matrix(Gfixed)
      if(nrow(Gfixed)!=nfac | ncol(Gfixed)!=nfac){stop("Input 'Gfixed' must have 'nfac' rows and 'nfac' columns")}
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
    
    # check 'Gstart' and 'Bstart' and 'Cstart' inputs
    if(!is.null(Gstart)){
      Gstart <- as.matrix(Gstart)
      if(nrow(Gstart)!=nfac | ncol(Gstart)!=nfac){stop("Input 'Gstart' must have 'nfac' rows and 'nfac' columns")}
    }
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
    
    # check 'Gstruc' and 'Bstruc' and 'Cstruc' inputs
    if(!is.null(Gstruc)){
      Gstruc <- as.matrix(Gstruc)
      if(nrow(Gstruc)!=nfac){stop("Input 'Gstruc' must have 'nfac' rows")}
      if(ncol(Gstruc)!=nfac){stop("Input 'Gstruc' must have 'nfac' columns")}
      if(any("logical"!=c(apply(Gstruc,1:2,class)))){stop("Input 'Gstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
    }
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
      stop("Input 'cl' must be cluster (created by makeCluster) when parallel=TRUE \n  See examples in documentation:  ?parafac2")
    }
    
    # parafac2 fitting
    if(lxdim==3L){
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac2_3way",data=X,xcx=xcx,
                              const=const,maxit=maxit,ctol=ctol,Gfixed=Gfixed,Bfixed=Bfixed,
                              Cfixed=Cfixed,Gstart=Gstart,Bstart=Bstart,Cstart=Cstart,
                              Gstruc=Gstruc,Bstruc=Bstruc,Cstruc=Cstruc)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac2_3way(data=X,nfac=nfac,xcx=xcx,const=const,maxit=maxit,
                                         ctol=ctol,Gfixed=Gfixed,Bfixed=Bfixed,Cfixed=Cfixed,
                                         Gstart=Gstart,Bstart=Bstart,Cstart=Cstart,
                                         Gstruc=Gstruc,Bstruc=Bstruc,Cstruc=Cstruc)
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
      if(!is.null(Dstruc)){
        Dstruc <- as.matrix(Dstruc)
        if(nrow(Dstruc)!=xdim[4]){stop("Input 'Dstruc' must have the same number of rows as dim(X)[4]")}
        if(ncol(Dstruc)!=nfac){stop("Input 'Dstruc' must have 'nfac' columns")}
        if(any("logical"!=c(apply(Dstruc,1:2,class)))){stop("Input 'Dstruc' must be a matrix with logical (TRUE/FALSE) entries.")}
      }
      if(parallel){
        nstartlist <- vector("list",nstart)
        nstartlist[1:nstart] <- nfac
        pfaclist <- parLapply(cl=cl,X=nstartlist,fun="parafac2_4way",data=X,xcx=xcx,
                              const=const,maxit=maxit,ctol=ctol,Gfixed=Gfixed,
                              Bfixed=Bfixed,Cfixed=Cfixed,Dfixed=Dfixed,
                              Gstart=Gstart,Bstart=Bstart,Cstart=Cstart,Dstart=Dstart,
                              Gstruc=Gstruc,Bstruc=Bstruc,Cstruc=Cstruc,Dstruc=Dstruc)
      } else {
        pfaclist <- vector("list",nstart)
        for(j in 1:nstart){
          pfaclist[[j]] <- parafac2_4way(data=X,nfac=nfac,xcx=xcx,const=const,maxit=maxit,ctol=ctol,
                                         Gfixed=Gfixed,Bfixed=Bfixed,Cfixed=Cfixed,Dfixed=Dfixed,
                                         Gstart=Gstart,Bstart=Bstart,Cstart=Cstart,Dstart=Dstart,
                                         Gstruc=Gstruc,Bstruc=Bstruc,Cstruc=Cstruc,Dstruc=Dstruc)
        }
      } # end if(parallel)
    } # end if(lxdim==3L) 
    
    # output results
    if(output[1]=="best"){
      widx <- which.max(sapply(pfaclist,function(x) x$Rsq))
      pfac <- pfaclist[[widx]]
      class(pfac) <- "parafac2"
      return(pfac)
    } else {
      pfaclist <- lapply(pfaclist, function(pfac) {
        class(pfac) <- "parafac2"
        pfac
        })
      return(pfaclist)
    }
    
  } # end parafac2.R