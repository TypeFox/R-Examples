qregcdf <- function(form,taumat=c(.10,.25,.50,.75,.90),hx=0,hy=0,nx=20,ny=100,targetx=0,targety=0,
    graph.target=FALSE,graph.yhat=FALSE, data=NULL) {

  xmat <- model.frame(form,data=data)
  namevect <- names(xmat)
  yname = namevect[1]
  xname = namevect[2]
  y <- xmat[,1]
  x <- xmat[,2]

  if (hx==0) {
    q25 = quantile(x,.25)
    q75 = quantile(x,.75)
    sx = sd(x)
    n = length(x)
    hx = 1.06*min(sx, (q75-q25)/1.349)*(n^(-.2))
  }

  if (hy==0) {
    q25 = quantile(y,.25)
    q75 = quantile(y,.75)
    sy = sd(y)
    n = length(y)
    hy = 1.06*min(sy, (q75-q25)/1.349)*(n^(-.2))
  }

  ntau = length(taumat)
  if (identical(targetx,0)) {targetx <- seq(min(x),max(x),length=nx) }
  nx = length(targetx)
  if (identical(targety,0)) {targety <- seq(min(y),max(y),length=ny) }
  ny = length(targety)
  fxy <- array(0,dim=c(nx,ny))

  for (i in 1:nx) {
  for (j in 1:ny) {
    kx <- dnorm((x-targetx[i])/hx)
    ky <- pnorm((targety[j]-y)/hy)
    fxy[i,j] = mean(kx*ky)/mean(kx)
  }
  }

  yhatmat.target <- array(0,dim=c(nx,ntau))

  for (i in 1:nx) {
    temp <- fxy[i,]
  for (j in 1:(ny-1)) {
  for (jtau in 1:ntau) {
      if (temp[j]<taumat[jtau] & temp[j+1]>taumat[jtau]) yhatmat.target[i,jtau] = targety[j]
  }
  }
  }


  ymin = min(yhatmat.target)
  ymax = max(yhatmat.target)
  if (graph.target==TRUE) {
    plot(targetx,yhatmat.target[,1],type="l",xlab=xname,ylab=yname,main="Predicted Values at Target Points",ylim=c(ymin,ymax))
    if (ntau>1) {
      for (j in seq(2,ntau)) {
        lines(targetx,yhatmat.target[,j])
      }
    }
  }

  n = length(y)
  yhatmat <- array(0,dim=c(n,ntau))
  for (j in seq(1,ntau)) {
    yhatmat[,j] <- smooth12(targetx,yhatmat.target[,j],x)
  }
  if (graph.yhat==TRUE) {
    plot(x,yhatmat[,1],type="l",xlab=xname,ylab=yname,main="Predicted Values",ylim=c(ymin,ymax))
    if (ntau>1) {
      for (j in seq(2,ntau)) {
        lines(x,yhatmat[,j])
      }
    }
  }
    
  out <- list(yhatmat,yhatmat.target,targetx,targety,taumat,hx,hy)
  names(out) <- c("yhat","yhat.target","targetx","targety","taumat","hx","hy")
  return(out)
}


