rfvar <-
function(ydata=NA, lags=2, xdata=NULL, const=TRUE, breaks=NULL){
  if(is.null(dim(ydata))) dim(ydata) <- c(length(ydata),1)
  T <- dim(ydata)[1]
  nvar <- dim(ydata)[2]
  if(const){
    xdata <- cbind(xdata,matrix(1,T,1))
  }
  nox <- identical(xdata,NULL)
  if(!nox){
    T2 <- dim(xdata)[1]
    nx <- dim(xdata)[2]
  }
  else{
    T2 <- T
    nx <- 0
    xdata <- matrix(0,T2,0)
  }
  #
  if(!identical(T2,T)){
    print('xdata and ydata must be same length.')
    return()
  }
  if(identical(breaks,NULL))
    nbreaks <- 0
  else
    nbreaks<-length(breaks)
  breaks <- c(0,breaks,T)
  if(any(breaks[2:length(breaks)] <= breaks[1:(length(breaks)-1)]))
    stop("list of breaks must be in strictly increasing order\n")
  if(breaks[2]>lags)
    smpl <- (lags+1):breaks[2]
  else
    smpl <- NULL
  if(nbreaks>0){
    for (nb in 2:(nbreaks+1))
      smpl <- c(smpl,(breaks[nb]+lags+1):breaks[nb+1])
  }
  Tsmpl <- length(smpl)
  X <- array(0,dim=c(Tsmpl,nvar,lags))
  for(ix in seq(along=smpl))
    X[ix,,] <- t(ydata[smpl[ix]-(1:lags),,drop=FALSE])
  dim(X) <- c(Tsmpl,nvar*lags)
  X <- cbind(X, xdata[smpl,,drop=FALSE]) # rhs data for model
  y <- ydata[smpl,,drop=FALSE] # lhs data for model
  vldvr <- svd(X)
  di <- 1./vldvr$d
  dfx <- sum(vldvr$d > 100*.Machine$double.eps)
  di <- di[1:dfx]
  vldvr$u <- vldvr$u[, 1:dfx]
  vldvr$v <- vldvr$v[, 1:dfx]
  snglty <- dim(X)[2] - dfx
  B <- vldvr$v %*% (di * (t(vldvr$u) %*% y))
  u <-  y-X %*% B
  if (!is.null(tsp(ydata))) u <- ts(u, start=start(ydata)+c(0,lags),frequency=frequency(ydata))
  nX <- dim(X)[2]
  xx <-  di * t(vldvr$v)
  xx <-  crossprod(xx)
  By <-  t(B) # fix this only plus 1 for const
  # dim(By) <-  c(nvar,lags,nvar)       # variables, lags, equations
  # By <-  aperm(By,c(3,1,2)) #equations, variables, lags to match impulsdt.m
  if(!is.null(dimnames(ydata)[2]))
  {
    ynames <- dimnames(ydata)[[2]]
  }else
  {
    ynames <- rep("",times=nvar)
  }
  if(!nox)
  {
    if(!is.null(dimnames(xdata)[2]))
    {
      xnames <- dimnames(xdata)[[2]]
      xxnames <- c(paste(rep(ynames,each=lags),1:lags,sep=""), xnames)
      # dimnames(xx) <- list(xxnames,xxnames)
      dimnames(By) <- list(ynames, c(paste(ynames,rep(1:lags, each=length(ynames)), sep="")))
    }else
    {
      xnames <- rep("",times=nx)
      xxnames <- c(paste(rep(ynames,each=lags),1:lags,sep=""))
      # dimnames(xx) <- list(xxnames,xxnames)
      dimnames(By) <- list(ynames, c(paste(ynames,rep(1:lags, each=length(ynames)), sep=""),"const"))
    }
  }
  if (nox)
    Bx <-  NULL
  else
  {
    Bx <-  matrix(B[nvar*lags+(1:nx),],dim(B)[2],nx)
    dimnames(Bx) <- list(ynames,xnames)
  }
  return(list(By=By, Bx=Bx, u=u, xx = xx, singular=snglty, X=X))
}
