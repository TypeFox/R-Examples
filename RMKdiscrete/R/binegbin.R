dbinegbin <- function(y, nu, p, log=FALSE, add.carefully=FALSE){ #nu & p are vectors of order 3; [0,1,2]  
  log <- log[1];  add.carefully <- add.carefully[1]
  if(length(dim(y))<2){y <- matrix(y,ncol=2,byrow=T)}
  if(length(dim(nu))<2){nu <- matrix(nu,ncol=3,byrow=T)}
  if(length(dim(p))<2){p <- matrix(p,ncol=3,byrow=T)}
  #print(c(nrow(x),nrow(nu),nrow(p)))
  numrows <- max( c( nrow(y), nrow(nu), nrow(p) ) )
  inmtx <- cbind( matrix(rep(t(y),length=2*numrows),nrow=numrows,ncol=2,byrow=T), 
                  matrix(rep(t(nu),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
                  matrix(rep(t(p),length=3*numrows),nrow=numrows,ncol=3,byrow=T)
  )
  nout <- nrow(inmtx)
  out <- rep(0,nout)
  inok <- rep(TRUE,nout)
  if(all( rowSums(!is.finite(inmtx))>0 )){return(rep(NA,nout))}
  if(any(!is.finite(inmtx))){
    out[rowSums(!is.finite(inmtx))>0] <- NA
    inok[rowSums(!is.finite(inmtx))>0] <- FALSE
    out[inok] <- dbinegbin(y=inmtx[inok,1:2],nu=inmtx[inok,3:5],p=inmtx[inok,6:8],log=log,
                        add.carefully=add.carefully)
    return(out)
  }
  if( any(round(inmtx[,1:2])!=inmtx[,1:2]) | any(inmtx[,1:2]<0) ){
    out[ .rowSumsWrapper(round(inmtx[,1:2])!=inmtx[,1:2])>0 | .rowSumsWrapper(inmtx[,1:2]<0)>0 ] <- 
      ifelse(log==T,-Inf,0)
    inok[ .rowSumsWrapper(round(inmtx[,1:2])!=inmtx[,1:2])>0 | .rowSumsWrapper(inmtx[,1:2]<0)>0 ] <- FALSE
  }
  if( any(inmtx[,3:5]<0) | any(inmtx[,6:8]>1) | any(inmtx[,6:8]<=0) ){
    warning("NaNs produced")
    out[ (.rowSumsWrapper(inmtx[,3:5]<0)>0) | (.rowSumsWrapper(abs(inmtx[,6:8])>1)>0) ] <- NaN
    inok[ (.rowSumsWrapper(inmtx[,3:5]<0)>0) | (.rowSumsWrapper(abs(inmtx[,6:8])>1)>0) ] <- FALSE
  }
  if(sum(!inok)<nout){
    #void call_dbinegbin(double *x, double *y, double *nu0, double *nu1, double *nu2, double *p0, double *p1, 
    #                    double *p2, int *give_log, int *add_carefully, int *Cnout, double *Cout){
    out[inok] <- .C( "call_dbinegbin",as.numeric(inmtx[inok,1]),as.numeric(inmtx[inok,2]),
               as.numeric(inmtx[inok,3]),as.numeric(inmtx[inok,4]),as.numeric(inmtx[inok,5]),
               as.numeric(inmtx[inok,6]),as.numeric(inmtx[inok,7]),as.numeric(inmtx[inok,8]),
               as.integer(log),as.integer(add.carefully),as.integer(sum(inok)),as.numeric(rep(0,sum(inok))) )[[12]]
  }
  return(out)
}


binegbin.logMV <- function(nu,p,const.add=1,tol=1e-14,add.carefully=FALSE){
  tol <- tol[1]; add.carefully <- add.carefully[1]
  if(length(dim(nu))<2){nu <- matrix(nu,ncol=3,byrow=T)}
  if(length(dim(p))<2){p <- matrix(p,ncol=3,byrow=T)}
  numrows <- max( c( nrow(nu), nrow(p), length(const.add), length(tol) ) )
  inmtx <- cbind( 
    matrix(rep(t(nu),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
    matrix(rep(t(p),length=3*numrows),nrow=numrows,ncol=3,byrow=T),
    rep(const.add,length=numrows),rep(tol,length=numrows)
  )
  nout <- nrow(inmtx)
  out <- matrix(NA,nrow=nout,ncol=5,dimnames=list(NULL,c("EY1","EY2","VY1","VY2","COV")))
  inok <- rep(TRUE,nout)
  if(all( rowSums(is.na(inmtx))>0 )){return(out)}
  if(any(is.na(inmtx))){
    inok[rowSums(is.na(inmtx))>0] <- FALSE
    out[inok,] <- binegbin.logMV(nu=inmtx[inok,1:3],p=inmtx[inok,4:6],const.add=inmtx[inok,7],tol=inmtx[inok,8])
    return(out)
  }
  if(any(inmtx[,7]<=0)){stop("argument 'const.add' must be greater than zero")}
  if(any(inmtx[,8]<=0)){stop("argument 'tol' must be greater than zero")}
  for(i in 1:nrow(out)){
    if(any(inmtx[i,1:3]<0) | any(inmtx[i,4:6]<=0) | any(inmtx[i,4:6]>1) | any(!is.finite(inmtx[i,1:6]))){
      out[i,] <- NaN
      next
    }
    else{
#       void call_binegbin_logMV(double *nu0, double *nu1, double *nu2, 
#                                double *p0, double *p1, double *p2,
#                                double *const_add, double *tol, int *add_carefully,
#                                double *EX, double *EY, double *EX2, double *EY2, double *EXY){
      out.temp <- .C("call_binegbin_logMV",as.double(inmtx[i,1]),as.double(inmtx[i,2]),as.double(inmtx[i,3]),
                     as.double(inmtx[i,4]), as.double(inmtx[i,5]), as.double(inmtx[i,6]),
                     as.double(inmtx[i,7]),as.double(inmtx[i,8]), as.integer(add.carefully),
                     as.double(0),as.double(0),as.double(0),as.double(0),as.double(0))
      out[i,] <- c(out.temp[[10]],
                   out.temp[[11]],
                   out.temp[[12]]-(out.temp[[10]]^2), #EX2-EX^2,
                   out.temp[[13]]-(out.temp[[11]]^2), #EY2-EY^2,
                   out.temp[[14]]-(out.temp[[10]]*out.temp[[11]])) #EXY-(EX*EY))
    }
  }
  return(out)
}

rbinegbin <- function(n, nu, p){ #nu & p: [0,1,2]  
  n <- as.integer(n[1])
  if(length(dim(nu))<2){nu <- matrix(nu,nrow=n,ncol=3,byrow=T)}
  else{nu <- matrix(t(nu),nrow=n,ncol=3,byrow=T)}
  if(length(dim(p))<2){p <- matrix(p,nrow=n,ncol=3,byrow=T)}
  else{p <- matrix(t(p),nrow=n,ncol=3,byrow=T)}
  u <- rnbinom(n=n,size=nu[,1],prob=p[,1])
  xp <- rnbinom(n=n,size=nu[,2],prob=p[,2])
  yp <- rnbinom(n=n,size=nu[,3],prob=p[,3])
  out <- cbind(u,u) + cbind(xp,yp)
  dimnames(out) <- NULL
  return(out)
}