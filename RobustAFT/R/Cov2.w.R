"Cov2.w" <-
function(X,y,l,u,theta,sigma,opt=c("integrals","averages"),xk=1.717817) {
n <- nrow(X); np <- ncol(X); n2 <- n*(n-np); xk <- 1.717817
AV0 <- matrix(0,ncol=np+1,nrow=np+1)
AV1 <- matrix(0,ncol=np+1,nrow=np+1)
rs    <- as.vector(y-X%*%as.matrix(theta))/sigma
wi    <- ((l<rs)&(rs<u))*1
sumwi <- sum(wi)
#XtX  <- (t(X)%*%diag(wi)%*%X)/sumwi
#xbar <- apply(wi*X,2,mean)*(n/sumwi)
  wiX <- wi*X 
  XtX <- (t(X) %*% wiX)/sumwi
 xbar <- apply(wiX, 2, mean) * (n/sumwi)  

#alfa <- Alpha.w(l,u); beta <- Beta.w(l,u)
l     <- max(l,-25); u <- min(u,4)
if (opt=="averages")  {invM0 <- invM2.w(l,u,theta,sigma,rs,wi,XtX,xbar,estim="SA")$Minv 
                       invM1 <- invM2.w(l,u,theta,sigma,rs,wi,XtX,xbar,estim="TMLA")$Minv}
if (opt=="integrals") {invM0 <- invM2.w(l,u,theta,sigma,rs,wi,XtX,xbar,estim="SI")$Minv 
                       invM1 <- invM2.w(l,u,theta,sigma,rs,wi,XtX,xbar,estim="TMLI")$Minv}
ncov  <- max(np*(np+1)/2,2)
avts0 <- matrix(double(1),np+1,np+1)
avts  <- matrix(double(1),np+1,np+1)
 f.res <- .Fortran("av_tmlwf",X=as.double(X),y=as.double(y),n=as.integer(n),
   np=as.integer(np),ncov=as.integer(ncov),l=as.double(l),u=as.double(u),
   xk=as.double(xk),theta=as.double(theta),sigma=as.double(sigma),
   invm0=as.double(invM0),invm1=as.double(invM1),avts0=as.double(avts0),
   avts=as.double(avts),xbar=as.double(xbar),XtX=as.double(XtX),
   sa=double(ncov),sc1=double(ncov),x0=double(np),its0=double(np+1),
   its=double(np+1))
 AV.TS0 <- matrix(f.res$avts0,nrow=np+1,ncol=np+1)
 AV.TS  <- matrix(f.res$avts,nrow=np+1,ncol=np+1)
list(CV0=AV.TS0, CV1=AV.TS, XtX=XtX, xbar=xbar)}

