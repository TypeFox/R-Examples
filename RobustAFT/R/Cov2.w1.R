"Cov2.w1" <-
function(y,l,u,theta,sigma,opt=c("integrals","averages")) {
n <- length(y); n2 <- n*(n-1); xk <- 1.717817
AV0 <- matrix(0,ncol=2,nrow=2)
AV1 <- matrix(0,ncol=2,nrow=2)
rs    <- (y-theta)/sigma
wi    <- ((l<rs)&(rs<u))*1
X     <- matrix(1,nrow=n,ncol=1)
XtX   <- as.matrix(1)
xbar  <- 1
#alfa <- Alpha.w(l,u); beta <- Beta.w(l,u)
l     <- max(l,-25); u <- min(u,4)
if (opt=="averages")  {invM0 <- invM2.w1(l,u,theta,sigma,rs,wi,estim="SA")$Minv 
                       invM1 <- invM2.w1(l,u,theta,sigma,rs,wi,estim="TMLA")$Minv}
if (opt=="integrals") {invM0 <- invM2.w1(l,u,theta,sigma,rs,wi,estim="SI")$Minv 
                       invM1 <- invM2.w1(l,u,theta,sigma,rs,wi,estim="TMLI")$Minv}
ncov  <- 2
avts0 <- matrix(double(1),2,2)
avts  <- matrix(double(1),2,2)
 f.res <- .Fortran("av_tmlwf",X=as.double(X),y=as.double(y),n=as.integer(n),
   np=as.integer(1),ncov=as.integer(ncov),l=as.double(l),u=as.double(u),
   xk=as.double(xk),theta=as.double(theta),sigma=as.double(sigma),
   invm0=as.double(invM0),invm1=as.double(invM1),avts0=as.double(avts0),
   avts=as.double(avts),xbar=as.double(xbar),XtX=as.double(XtX),
   sa=double(ncov),sc1=double(ncov),x0=double(1),its0=double(2),its=double(2))
 AV.TS0 <- matrix(f.res$avts0,nrow=2,ncol=2)
 AV.TS  <- matrix(f.res$avts,nrow=2,ncol=2)
list(CV0=AV.TS0, CV1=AV.TS)}

