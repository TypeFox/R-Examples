"Cov2.n1" <-
function(y,u,theta,sigma,opt=c("integrals","averages")){
n <- length(y); np <- 1; n2 <- n*(n-1); xk <- 1.5477
AV0 <- matrix(0,ncol=2,nrow=2)
AV1 <- matrix(0,ncol=2,nrow=2)
rs    <- (y-theta)/sigma
wi    <- (abs(rs)<u)*1
X     <- matrix(1,nrow=n,ncol=1)
XtX   <- as.matrix(1)
xbar  <- 1
#alfa <- Alpha.n(u); beta <- Beta.n(u)
u     <- min(u,10)
if (opt=="averages")  {invM0 <- invM2.n1(u,theta,sigma,rs,wi,estim="SA")$Minv 
                       invM1 <- invM2.n1(u,theta,sigma,rs,wi,estim="TMLA")$Minv}
if (opt=="integrals") {invM0 <- invM2.n1(u,theta,sigma,rs,wi,estim="SI")$Minv 
                       invM1 <- invM2.n1(u,theta,sigma,rs,wi,estim="TMLI")$Minv}
ncov  <- 2
avts0 <- matrix(double(1),2,2)
avts  <- matrix(double(1),2,2)
f.res <- .Fortran("av_tmlnf",X=as.double(X),y=as.double(y),n=as.integer(n),
   np=as.integer(np),ncov=as.integer(ncov),u=as.double(u),k0=as.double(xk),
   theta=as.double(theta),sigma=as.double(sigma),invm0=as.double(invM0),
   invm1=as.double(invM1),avts0=as.double(avts0),avts=as.double(avts),
   xbar=as.double(xbar),XtX=as.double(XtX),sa=double(ncov),sc1=double(ncov),
   x0=double(np),its0=double(2),its=double(2))
 AV.TS0 <- matrix(f.res$avts0,nrow=2,ncol=2)
 AV.TS  <- matrix(f.res$avts,nrow=2,ncol=2)
list(CV0=AV.TS0, CV1=AV.TS)}

