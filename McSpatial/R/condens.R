condens <- function(form,window=.7,bandwidth=0,kern="tcub",mingrid.x=NULL,maxgrid.x=NULL,mingrid.y=NULL,maxgrid.y=NULL,ngrid=50,
  xlab="x",ylab="y",zlab="fxy/fx",
  contour=TRUE,level=TRUE,wire=TRUE,dens=TRUE,targetx.dens=NULL,quantile.dens=c(.10,.25,.50,.75,.90),data=NULL) {

  library(lattice)
  mat <- model.frame(form,data=data)
  y <- mat[,1]
  x <- mat[,2]
  n = length(x)

  if (kern=="rect")  { wgt <- function(psi) {ifelse(abs(psi)>=0,.5,0) } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { .75*(1-psi^2) } }
  if (kern=="bisq")  { wgt <- function(psi) { (15/16)*((1-psi^2)^2) } }
  if (kern=="tcub")  { wgt <- function(psi) { (70/81)*((1 - abs(psi)^3)^3) } }
  if (kern=="trwt")  { wgt <- function(psi) { (35/32)*((1 - psi^2)^3) } }
  if (kern=="gauss") { wgt <- function(psi) { dnorm(psi) } }

  if (window>0){
    targetx <- maketarget(~x,window=window,bandwidth=0,kern=kern)$target
    targety <- maketarget(~y,window=window,bandwidth=0,kern=kern)$target
  }
  if (bandwidth>0){
    targetx <- maketarget(~x,window=0,bandwidth=bandwidth*sd(x),kern=kern)$target
    targety <- maketarget(~y,window=0,bandwidth=bandwidth*sd(y),kern=kern)$target
  }
  targetx <- sort(targetx)
  targety <- sort(targety)

  targetxy <- expand.grid(targetx,targety)
  nx = length(targetx)
  ny = length(targety)
  nxy = nrow(targetxy)
  fx.target  <- array(0,dim=nx)
  fy.target  <- array(0,dim=ny)
  fxy.target <- array(0,dim=nxy)
  
  for (i in seq(1:nx)) {
    dist <- abs(x-targetx[i])
    maxd = ifelse(bandwidth==0,quantile(dist,window),bandwidth*sd(x))
    k <- wgt(dist/maxd)
    k <- ifelse(kern!="gauss"&dist>maxd,0,k)
    fx.target[i] = mean(k)/maxd
  }
  for (i in seq(1:ny)) {
    dist <- abs(y-targety[i])
    maxd = ifelse(bandwidth==0,quantile(dist,window),bandwidth*sd(y))
    k <- wgt(dist/maxd)
    k <- ifelse(kern!="gauss"&dist>maxd,0,k)
    fy.target[i] = mean(k)/maxd
  }
  for (i in seq(1:nxy)) {
    dist <- abs(x-targetxy[i,1])
    maxd = ifelse(bandwidth==0,quantile(dist,window),bandwidth*sd(x))
    kx <- wgt(dist/maxd)
    kx <- ifelse(kern!="gauss"&dist>maxd,0,kx)
    kx <- kx/maxd

    dist <- abs(y-targetxy[i,2])
    maxd = ifelse(bandwidth==0,quantile(dist,window),bandwidth*sd(y))
    ky <- wgt(dist/maxd)
    ky <- ifelse(kern!="gauss"&dist>maxd,0,ky)
    ky <- ky/maxd

    fxy.target[i] = mean(kx*ky)
  }

  fx <- smooth12(targetx,fx.target,x)
  fy <- smooth12(targety,fy.target,y)
  fxy <- smooth12(targetxy,fxy.target,cbind(x,y))


  if (identical(mingrid.x,NULL)) { mingrid.x = min(x) }
  if (identical(maxgrid.x,NULL)) { maxgrid.x = max(x) }
  if (identical(mingrid.y,NULL)) { mingrid.y = min(y) }
  if (identical(maxgrid.y,NULL)) { maxgrid.y = max(y) }

  grid.x <- seq(mingrid.x,maxgrid.x,length=ngrid) 
  grid.y <- seq(mingrid.y,maxgrid.y,length=ngrid)
  xy <- expand.grid(grid.x,grid.y)
  grid.x <- xy[,1]
  grid.y <- xy[,2]
  grid.fxy <- smooth12(targetxy,fxy.target,xy)
  grid.fx <- smooth12(targetx,fx.target,grid.x)
  grid.fxy <- grid.fxy/grid.fx
  gridmat <- cbind(grid.x,grid.y,grid.fxy)

  if (contour==TRUE) {print(contourplot(grid.fxy~grid.x*grid.y,xlab=xlab,ylab=ylab))}
  if (level==TRUE)     {print(levelplot(grid.fxy~grid.x*grid.y,xlab=xlab,ylab=ylab))}
  if (wire==TRUE)      {print(wireframe(grid.fxy~grid.x*grid.y,xlab=xlab,ylab=ylab,zlab=zlab))}

  if (identical(targetx.dens,NULL)) {targetx.dens <- quantile(x,quantile.dens) }
  if (length(targetx.dens)>5) {
    cat("Target x for density graph > 5;  will only use first 5","\n")
    targetx.dens <- targetx.dens[1:5]
  }
  nq <- length(targetx.dens)
  densmat <- array(0,dim=c(n,nq))
  for (j in seq(1:nq)) {
    dist <- abs(x-targetx.dens[j])
    maxd = ifelse(bandwidth==0,quantile(dist,window),bandwidth*sd(x))
    kx <- wgt(dist/maxd)
    kx <- ifelse(kern!="gauss"&dist>maxd,0,kx)/maxd
    fx1 = mean(kx)

    for (i in seq(1:ny)) {
      dist <- abs(y-targety[i])
      maxd = ifelse(bandwidth==0,quantile(dist,window),bandwidth*sd(y))
      ky <- wgt(dist/maxd)
      ky <- ifelse(kern!="gauss"&dist>maxd,0,ky)/maxd
      fy.target[i] = mean(kx*ky)/fx1
    }
    densmat[,j] <- smooth12(targety,fy.target,y)

  }
    
  if (dens==TRUE) {
    colmat <- c("black","blue","red","green","orange")
    o <- order(y)
    plot(y[o],densmat[o,1],xlab=ylab,ylab=zlab,type="l",ylim=c(min(densmat,na.rm=TRUE),max(densmat,na.rm=TRUE)))
    for (j in seq(2,nq)) {
      lines(y[o],densmat[o,j],col=colmat[j])
    }
    legend("topright",as.character(targetx.dens),col=colmat[1:nq],lwd=1) 
  }

  out <- list(fx,fy,fxy,gridmat,densmat)
  names(out) <- c("fx","fy","fxy","gridmat","densmat")
  return(out)
}
