.rowSumsWrapper <- function(x,na.rm=F){
  if(is.null(dim(x))){return(sum(x,na.rm=na.rm))}
  else{return(rowSums(x,na.rm=na.rm))}
}

dbiLGP <- function(y, theta, lambda, nc=NULL, log=FALSE, add.carefully=FALSE){ #theta & lambda are vectors of order 3; [0,1,2]  
  log <- log[1];  add.carefully <- add.carefully[1]
  if(length(dim(y))<2){y <- matrix(y,ncol=2,byrow=T)}
  if(length(dim(theta))<2){theta <- matrix(theta,ncol=3,byrow=T)}
  if(length(dim(lambda))<2){lambda <- matrix(lambda,ncol=3,byrow=T)}
  if(is.null(nc)){
    nc <- cbind(LGP.get.nc(theta=theta[,1],lambda=lambda[,1]),LGP.get.nc(theta=theta[,2],lambda=lambda[,2]),
                LGP.get.nc(theta=theta[,3],lambda=lambda[,3]))
  }
  else{ if(length(dim(nc))<2){nc <- matrix(nc,ncol=3,byrow=T)} }
  numrows <- max( c( nrow(y), nrow(theta), nrow(lambda), nrow(nc) ) )
  inmtx <- cbind( matrix(rep(t(y),length=2*numrows),nrow=numrows,ncol=2,byrow=T), 
                  matrix(rep(t(theta),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
                  matrix(rep(t(lambda),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
                  matrix(rep(t(nc),length=3*numrows),nrow=numrows,ncol=3,byrow=T)
  )
  nout <- nrow(inmtx)
  out <- rep(0,nout)
  inok <- rep(TRUE,nout)
  if(all( rowSums(!is.finite(inmtx))>0 )){return(rep(NA,nout))}
  if(any(!is.finite(inmtx))){
    out[rowSums(!is.finite(inmtx))>0] <- NA
    inok[rowSums(!is.finite(inmtx))>0] <- FALSE
    out[inok] <- dbiLGP(y=inmtx[inok,1:2],theta=inmtx[inok,3:5],lambda=inmtx[inok,6:8],nc=inmtx[inok,9:11],log=log,
                        add.carefully=add.carefully)
    return(out)
  }
  if( any(round(inmtx[,1:2])!=inmtx[,1:2]) | any(inmtx[,1:2]<0) ){
    out[ .rowSumsWrapper(round(inmtx[,1:2])!=inmtx[,1:2])>0 | .rowSumsWrapper(inmtx[,1:2]<0)>0 ] <- 
      ifelse(log==T,-Inf,0)
    inok[ .rowSumsWrapper(round(inmtx[,1:2])!=inmtx[,1:2])>0 | .rowSumsWrapper(inmtx[,1:2]<0)>0 ] <- FALSE
  }
  if( any(inmtx[,3:5]<0) | any(abs(inmtx[,6:8])>1) ){
    warning("NaNs produced")
    out[ (.rowSumsWrapper(inmtx[,3:5]<0)>0) | (.rowSumsWrapper(abs(inmtx[,6:8])>1)>0) ] <- NaN
    inok[ (.rowSumsWrapper(inmtx[,3:5]<0)>0) | (.rowSumsWrapper(abs(inmtx[,6:8])>1)>0) ] <- FALSE
  }
  if(sum(!inok)<nout){
    #void call_dbiLGP(int *x, int *y, double *theta0, double *theta1, double *theta2, double *lambda0, double *lambda1,
    #double *lambda2, double *nc0, double *nc1, double *nc2, int *give_log, int *add_carefully int *Cnout, double *Cout){
    out[inok] <- .C("call_dbiLGP",as.numeric(inmtx[inok,1]),as.numeric(inmtx[inok,2]),as.numeric(inmtx[inok,3]),
                    as.numeric(inmtx[inok,4]),as.numeric(inmtx[inok,5]),as.numeric(inmtx[inok,6]),
                    as.numeric(inmtx[inok,7]),as.numeric(inmtx[inok,8]),as.numeric(inmtx[inok,9]),
                    as.numeric(inmtx[inok,10]),as.numeric(inmtx[inok,11]),as.integer(log),as.integer(add.carefully),
                    as.integer(sum(inok)),as.double(rep(0,sum(inok))))[[15]]
  }
  return(out)
}

biLGP.logMV <- function(theta,lambda,nc=NULL,const.add=1,tol=1e-14,add.carefully=FALSE){
  tol <- tol[1]; add.carefully <- add.carefully[1]
  if(length(dim(theta))<2){theta <- matrix(theta,ncol=3,byrow=T)}
  if(length(dim(lambda))<2){lambda <- matrix(lambda,ncol=3,byrow=T)}
  if(is.null(nc)){
    nc <- cbind(LGP.get.nc(theta=theta[,1],lambda=lambda[,1]),LGP.get.nc(theta=theta[,2],lambda=lambda[,2]),
                LGP.get.nc(theta=theta[,3],lambda=lambda[,3]))
  }
  else{ if(length(dim(nc))<2){nc <- matrix(nc,ncol=3,byrow=T)} }
  numrows <- max( c( nrow(theta), nrow(lambda), length(const.add), nrow(nc) ) )
  inmtx <- cbind( 
    matrix(rep(t(theta),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
    matrix(rep(t(lambda),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
    matrix(rep(t(nc),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
    rep(const.add,length=numrows)
  )
  nout <- nrow(inmtx)
  out <- matrix(NA,nrow=nout,ncol=5,dimnames=list(NULL,c("EY1","EY2","VY1","VY2","COV")))
  inok <- rep(TRUE,nout)
  if(all( rowSums(!is.finite(inmtx))>0 )){return(out)}
  if(any(!is.finite(inmtx))){
    inok[rowSums(!is.finite(inmtx))>0] <- FALSE
    out[inok,] <- biLGP.logMV(theta=inmtx[inok,1:3],lambda=inmtx[inok,4:6],nc=inmtx[inok,7:9],
                              const.add=inmtx[inok,10],tol=tol,add.carefully=add.carefully)
    return(out)
  }
  if(any(inmtx[,10]<=0)){stop("argument 'const.add' must be greater than zero")}
  if(tol<=0){stop("argument 'tol' must be greater than zero")}
  if( any(inmtx[,1:3]<0) | any(abs(inmtx[,4:6])>1) ){
    warning("NaNs produced")
    out[ (.rowSumsWrapper(inmtx[,1:3]<0)>0) | (.rowSumsWrapper(abs(inmtx[,4:6])>1)>0) ] <- NaN
    inok[ (.rowSumsWrapper(inmtx[,1:3]<0)>0) | (.rowSumsWrapper(abs(inmtx[,4:6])>1)>0) ] <- FALSE
    out[inok,] <- biLGP.logMV(theta=inmtx[inok,1:3],lambda=inmtx[inok,4:6],nc=inmtx[inok,7:9],
                              const.add=inmtx[inok,10],tol=tol,add.carefully=add.carefully)
    return(out)
  }
  for(i in 1:nrow(out)){
#     void call_biLGP_logMV(double *theta0, double *theta1, double *theta2, 
#                           double *lambda0, double *lambda1, double *lambda2,
#                           double *nc0, double *nc1, double *nc2,
#                           double *const_add, double *tol, int *add_carefully,
#                           double *EX, double *EY, double *EX2, double *EY2, double *EXY)
    out.temp <- .C("call_biLGP_logMV",as.double(inmtx[i,1]),as.double(inmtx[i,2]),as.double(inmtx[i,3]),
                   as.double(inmtx[i,4]), as.double(inmtx[i,5]), as.double(inmtx[i,6]),
                   as.double(inmtx[i,7]),as.double(inmtx[i,8]),as.double(inmtx[i,9]),
                   as.double(inmtx[i,10]),as.double(tol),as.integer(add.carefully),
                   as.double(0),as.double(0),as.double(0),as.double(0),as.double(0))
    out[i,] <- c(out.temp[[13]],
                 out.temp[[14]],
                 out.temp[[15]]-(out.temp[[13]]^2), #EX2-EX^2,
                 out.temp[[16]]-(out.temp[[14]]^2), #EY2-EY^2,
                 out.temp[[17]]-(out.temp[[13]]*out.temp[[14]])) #EXY-(EX*EY))
  }
  return(out)
}

rbiLGP <- function(n, theta, lambda){ #theta & lambda: [0,1,2]  
  n <- as.integer(n[1])
  if(length(dim(theta))<2){theta <- matrix(theta,nrow=n,ncol=3,byrow=T)}
  else{theta <- matrix(t(theta),nrow=n,ncol=3,byrow=T)}
  if(length(dim(lambda))<2){lambda <- matrix(lambda,nrow=n,ncol=3,byrow=T)}
  else{lambda <- matrix(t(lambda),nrow=n,ncol=3,byrow=T)}
  u <- rLGP(n=n,theta=theta[,1],lambda=lambda[,1])
  xp <- rLGP(n=n,theta=theta[,2],lambda=lambda[,2])
  yp <- rLGP(n=n,theta=theta[,3],lambda=lambda[,3])
  out <- cbind(u,u) + cbind(xp,yp)
  dimnames(out) <- NULL
  return(out)
}
