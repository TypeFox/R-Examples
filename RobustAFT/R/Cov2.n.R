"Cov2.n" <-
function(X,y,u,theta,sigma,opt=c("integrals","averages"),xk=1.5477){
n <- nrow(X); np <- ncol(X); n2 <- n*(n-np)
AV0 <- matrix(0,ncol=np+1,nrow=np+1)
AV1 <- matrix(0,ncol=np+1,nrow=np+1)
rs    <- as.vector(y-X%*%theta)/sigma
wi    <- (abs(rs)<u)*1
sumwi <- sum(wi)
#XtX  <- (t(X)%*%diag(wi)%*%X)/sumwi
#xbar <- apply(wi*X,2,mean)*(n/sumwi) 
  wiX <- wi*X 
  XtX <- (t(X) %*% wiX)/sumwi 
 xbar <- apply(wiX, 2, mean) * (n/sumwi)  

#alfa <- Alpha.n(u); beta <- Beta.n(u)
u     <- min(u,10)
if (opt=="averages")  {invM0 <- invM2.n(u,theta,sigma,rs,wi,XtX,xbar,estim="SA")$Minv 
                       invM1 <- invM2.n(u,theta,sigma,rs,wi,XtX,xbar,estim="TMLA")$Minv}
if (opt=="integrals") {invM0 <- invM2.n(u,theta,sigma,rs,wi,XtX,xbar,estim="SI")$Minv 
                       invM1 <- invM2.n(u,theta,sigma,rs,wi,XtX,xbar,estim="TMLI")$Minv}
ncov  <- max(np*(np+1)/2,2)
avts0 <- matrix(double(1),np+1,np+1)
avts  <- matrix(double(1),np+1,np+1)

 f.res <- .Fortran("av_tmlnf",X=as.double(X),y=as.double(y),n=as.integer(n),
   np=as.integer(np),ncov=as.integer(ncov),u=as.double(u),k0=as.double(xk),
   theta=as.double(theta),sigma=as.double(sigma),invm0=as.double(invM0),
   invm1=as.double(invM1),avts0=as.double(avts0),avts=as.double(avts),
   xbar=as.double(xbar),XtX=as.double(XtX),sa=double(ncov),sc1=double(ncov),
   x0=double(np),its0=double(np+1),its=double(np+1))
 AV.TS0 <- matrix(f.res$avts0,nrow=np+1,ncol=np+1)
 AV.TS  <- matrix(f.res$avts,nrow=np+1,ncol=np+1)
list(CV0=AV.TS0, CV1=AV.TS, XtX=XtX, xbar=xbar)}

