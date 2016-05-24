ngca <- function(data,L=c(1000,1000,1000),T=10,m=3,eps=1.5,npca=min(dim(x)[2],dim(x)[1]),filter.time="None",filter.space=FALSE,method="temporal",dg.trend=2,h.space=3,h.time=3,keepv=TRUE,delta=NULL){
#
#  NGCA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(all(class(data)=="fmridata")) {
   data <- fmri.detrend(data,dg.trend)
   x <- extract.data(data)
   mask <- data$mask>0
   fmriobj <- TRUE
# extracted data in x
#
#   we should also remove trends here
#
   if(filter.time %in% c("High","Both")){
# high pass filter in time if specified
      d <- dim(x)[4]
      x <- x[,,,-1]-x[,,,-d]
   }
   if(filter.time %in% c("Low","Both")){
# low pass filter in time if specified
      dx <- dim(x)
      cat("Start smoothing in time (Bandwidth=",h.time,")\n")
      x <- .Fortran("smtime",
                    as.double(x),
                    as.integer(dx[1]),
                    as.integer(dx[2]),
                    as.integer(dx[3]),
                    as.integer(dx[4]),
                    as.logical(mask),
                    as.double(h.time),
                    xnew=double(prod(dx)),
                    double(as.integer(2*h.time+1)),
                    as.integer(2*h.time+1),
                    DUP=TRUE,
                    PACKAGE="fmri")$xnew
     cat("Smoothing in time finished)\n")
     dim(x) <- dx
   }
   if(filter.space) {
# spatial smoothing if specified
      dx <- dim(x)
      cat("Start spatial smoothing (Bandwidth=",h.space,")\n")
      wghts <- data$weights
      x <- .Fortran("smspace",
                    as.double(x),
                    as.integer(dx[1]),
                    as.integer(dx[2]),
                    as.integer(dx[3]),
                    as.integer(dx[4]),
                    as.logical(mask),
                    as.double(h.space),
                    xnew=double(prod(dx)),
                    as.double(wghts),
                    double(prod(as.integer(2*h.space/wghts+1))),
                    as.integer(as.integer(2*h.space/wghts[1]+1)),
                    as.integer(as.integer(2*h.space/wghts[2]+1)),
                    as.integer(as.integer(2*h.space/wghts[3]+1)),
                    DUP=TRUE,
                    PACKAGE="fmri")$xnew
     cat("Spatial moothing finished)\n")
     dim(x) <- dx
   }
} else if (class(data)%in%c("matrix","array")){
   x <- data
   mask <- TRUE
   fmriobj <- FALSE
} else {
   warning("data has incompatible class argument")
   return(data)
}
#  
#  x - data matrix  (Nxd)
#
cat("Start Non-Gaussian Component Analysis\n")
set.seed(1)
xdim <- dim(x)
lxdim <- length(xdim)
d <- xdim[lxdim]
n <- nn <- prod(xdim[1:(lxdim-1)])
dim(x) <- c(n,d)
mask <- as.vector(mask)
if(length(mask)==1) mask <- rep(mask,n)
x <- x[mask,]
n <- sum(mask)
if(is.null(npca)||npca >= min(d,n)) npca <- min(d,n)
# sweep out the mean  
xmean <- switch(method,"spatial"=apply(x,1,mean),
                       "temporal"=apply(x,2,mean))
x <- switch(method,"spatial"=sweep(x,1,xmean),
            "temporal"=sweep(x,2,xmean))
if(method=="temporal") x <- sweep(x,1,apply(x,1,mean))
if(method=="spatial"){
x <- t(x)
d <- dim(x)[2]
n <- dim(x)[1]
# this means d>n
}
svdx <- svd(x,nu=npca,nv=npca)
# reduce to space of first npca components
nd <- switch(method,"spatial"=d-1,"temporal"=n-1)
y <- svdx$v[1:npca,]%*%t(svdx$u)*sqrt(nd)
sigmainvhalf <- svdx$v[1:npca,]%*%diag(svdx$d[1:npca]^(-1))%*%t(svdx$v[1:npca,])*sqrt(nd)
sigmahalf <- svdx$v[1:npca,]%*%diag(svdx$d[1:npca])%*%t(svdx$v[1:npca,])/sqrt(nd)
#
#
Lsum <- L[1]+L[2]+2*L[3]
s <- c(if(L[1]>0) seq(.5,5,length=L[1]), 
       if(L[2]>0) seq(0.05,5,length=L[2]), 
       if(L[3]>0) seq(0.05,4,length=L[3]),
       if(L[3]>0) seq(0.05,4,length=L[3]))
Lsum <- L[1]+L[2]+2*L[3]
ifun <- c(rep(1,L[1]),rep(2,L[2]),rep(3,L[3]),rep(4,L[3]))
#
#   now fast ICA
#
if(is.null(delta)) {
    delta <- npca*5e-4
}
cat("delta set to",delta,"\n")
omega <- matrix(rnorm(Lsum*npca),npca,Lsum)
omega <- sweep(omega,2,sqrt(apply(omega^2,2,sum)),"/")
fz <- .Fortran("fastica",
              as.double(t(y)),
              as.double(omega),
              as.integer(npca),
              as.integer(n),
              as.integer(Lsum),
              as.integer(ifun),
              as.integer(T),
              double(npca),
              v=double(npca*Lsum),
              normv=double(Lsum),
              as.double(s),
              double(npca),
              as.double(delta),
              DUP=TRUE,
              PACKAGE="fmri")[c("v","normv")]
v <- fz$v
normv <- fz$normv
plot(density(normv))
if(eps>quantile(normv,.95)) {
eps <- quantile(normv,.95)
cat("eps reset to",eps,"\n")
}
dim(v) <- c(npca,Lsum)
v0 <- v
v <- t(v[,normv>eps])
jhat <- prcomp(v)
if(method=="spatial") {
xhat <- sigmahalf[,1:npca]%*%jhat$rotation[,1:m]
ihat <- x[,1:npca]%*%xhat 
} else {
ihat <- sigmahalf[,1:npca]%*%jhat$rotation[,1:m]
xhat <- x%*%ihat
}
if(fmriobj){
z <- matrix(0,nn,m)
z[mask,]<-xhat
xhat <- array(z,c(xdim[1:3],m))
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat)
if(keepv) {
z$v <- v0
z$normv <- normv
z$svdx <- svdx
z$xmean <- xmean
z$npca <- npca
z$method <- method
z$sinvhalf <- sigmainvhalf
z$shalf <- sigmahalf
z$mask <- mask
z$xdim <- xdim
z$method <- z$method
z$npca <- npca
}
class(z) <- "fmringca"
} else {
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat,jhat=jhat,svdx=svdx)
if(keepv) {
z$v <- v0
z$normv <- normv
}
class(z) <- "ngca"
}
z
}
ngca.thresh <- function(ngcaobj,eps=1.5,m=3){
#
#  NGCA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(class(ngcaobj)!="fmringca"||is.null(ngcaobj$normv)) return(ngcaobj)
method <- ngcaobj$method
if(is.null(method)) method <- "temporal"
# just to be able to check with previous calculations
v <- ngcaobj$v
normv <- ngcaobj$normv
v <- t(v[,normv>eps])
jhat <- prcomp(v)
svdx <- ngcaobj$svdx
npca <- ngcaobj$npca
if(is.null(npca)) npca <- length(svdx$d)
# just to be able to check with previous calculations
x <- svdx$u%*%diag(svdx$d)%*%t(svdx$v)
if(method=="spatial") {
xhat <- ngcaobj$shalf[,1:npca]%*%jhat$rotation[,1:m]
ihat <- x%*%ihat
} else {
ihat <- ngcaobj$shalf[,1:npca]%*%jhat$rotation[,1:m]
xhat <- x%*%ihat
}
mask <- ngcaobj$mask
nn <- prod(ngcaobj$xdim[1:3])
z <- matrix(0,nn,m)
z[mask,]<-xhat
xhat <- array(z,c(ngcaobj$xdim[1:3],m))
z <- list(ihat=ihat,sdev=jhat$sdev[1:m],xhat=xhat,method=method,npca=npca)
class(z) <- "fmringca"
z
}

ngcafmri.plot <- function(ngcaobj,qeps=.9,m=1,center=FALSE,scale=FALSE){
#
#  NGCA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(class(ngcaobj)!="fmringca"||is.null(ngcaobj$normv)) return(ngcaobj)
method <- ngcaobj$method
if(is.null(method)) method <- "temporal"
# just to be able to check with previous calculations
v <- ngcaobj$v
normv <- ngcaobj$normv
eps <- quantile(normv,qeps)
cat("Using threshold of",eps,"\n")
v <- t(v[,normv>eps])
jhat <- prcomp(v,center=center,scale=scale)

svdx <- ngcaobj$svdx
npca <- ngcaobj$npca
if(is.null(npca)) npca <- length(svdx$d)
# just to be able to check with previous calculations
x <- svdx$u%*%diag(svdx$d)%*%t(svdx$v)
if(method=="spatial") {
xhat <- ngcaobj$shalf[,1:npca]%*%jhat$rotation[,1:m]
ihat <- x%*%ihat
} else {
ihat <- ngcaobj$shalf[,1:npca]%*%jhat$rotation[,1:m]
xhat <- x%*%ihat
}
mask <- ngcaobj$mask
xdim <- ngcaobj$xdim
nn <- prod(xdim[1:3])
chat <- numeric(nn)
chat[mask]<-x%*%ihat[,m]
dim(chat) <- xdim[1:3]
dd <- xdim[3]
nd1 <- as.integer(sqrt(dd))+1
nd2 <- as.integer((dd)/nd1)+1
par(mfrow=c(nd1,nd2),mar=c(.5,.5,.5,.1))
#
#  now compute correlations
#
for(i in 1:dd) image(chat[,,i],zlim=range(chat))
plot(ihat[,m],type="l")
invisible(NULL)
}

fmriica <- function(data,m=3,method="temporal",xind=NULL,yind=NULL,zind=NULL,tind=NULL,
filter.time="None",filter.space=FALSE,h.space=3,h.time=3,keepv=FALSE,...){
#
#  ICA algorithm  for fMRI
#  x should be either a fMRI object or a matrix 
#
if(!require(fastICA)) stop("Please install package fastICA from CRAN")
if(all(class(data)=="fmridata")) {
   x <- extract.data(data)
   mask <- data$mask>0
}  else {
   warning("data has incompatible class argument")
   return(data)
}
if(!is.null(xind)) mask[-xind,,] <- FALSE
if(!is.null(yind)) mask[,-yind,] <- FALSE
if(!is.null(zind)) mask[,,-zind] <- FALSE
   if(filter.time %in% c("Low","Both")){
      dx <- dim(x)
      cat("Start smoothing in time (Bandwidth=",h.time,")\n")
      x <- .Fortran("smtime",
                    as.double(x),
                    as.integer(dx[1]),
                    as.integer(dx[2]),
                    as.integer(dx[3]),
                    as.integer(dx[4]),
                    as.logical(mask),
                    as.double(h.time),
                    xnew=double(prod(dx)),
                    double(as.integer(2*h.time+1)),
                    as.integer(2*h.time+1),
                    DUP=TRUE,
                    PACKAGE="fmri")$xnew
     cat("Smoothing in time finished)\n")
     dim(x) <- dx
   }
   if(filter.space) {
      dx <- dim(x)
      cat("Start spatial smoothing (Bandwidth=",h.space,")\n")
      wghts <- data$weights
      x <- .Fortran("smspace",
                    as.double(x),
                    as.integer(dx[1]),
                    as.integer(dx[2]),
                    as.integer(dx[3]),
                    as.integer(dx[4]),
                    as.logical(mask),
                    as.double(h.space),
                    xnew=double(prod(dx)),
                    as.double(wghts),
                    double(prod(as.integer(2*h.space/wghts+1))),
                    as.integer(as.integer(2*h.space/wghts[1]+1)),
                    as.integer(as.integer(2*h.space/wghts[2]+1)),
                    as.integer(as.integer(2*h.space/wghts[3]+1)),
                    DUP=TRUE,
                    PACKAGE="fmri")$xnew
     cat("Spatial moothing finished)\n")
     dim(x) <- dx
   }
#
#  x - data matrix  (Nxd)
#
set.seed(1)
xdim <- dim(x)
lxdim <- length(xdim)
d <- dd <- xdim[lxdim]
n <- nn <- prod(xdim[1:(lxdim-1)])
if(is.null(tind)) tind <- (1:d)
dim(x) <- c(n,d)
mask <- as.vector(mask)
if(length(mask)==1) mask <- rep(mask,n)
x <- x[mask,tind]
n <- sum(mask)
if(method=="spatial"){
x <- t(x)
d <- dim(x)[2]
n <- dim(x)[1]
}
#
#   now fast ICA
#
ttt <- fastICA(x,m,method="C",...)
if(method=="spatial"){
z <- matrix(0,nn,m)
z[mask,]<-ttt$A
ihat <- array(z,c(xdim[1:3],m))
xhat <- matrix(0,m,dd)
xhat[,tind]<-ttt$S
} else {
z <- matrix(0,nn,m)
z[mask,]<-ttt$S
xhat <- array(z,c(xdim[1:3],m))
ihat <- matrix(0,m,dd)
ihat[,tind] <- ttt$A
}
z <- list(ihat=ihat,xhat=xhat)
class(z) <- "fmringca"
z
}

