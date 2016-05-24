edr <- function(x,y,m=2,rho0=1,h0=NULL,ch=exp(.5/max(4,(dim(x)[2]))),
    crhomin=1,cm=4,method="Penalized",basis="Quadratic",cw=NULL,graph=FALSE,
    show=1,trace=FALSE,fx=NULL,R=NULL){
args <- match.call()      
expand <- function(x) {
    d <- length(x)
    outer(x,x)[rep(1:d,d)>=rep(1:d,rep(d,d))]
    }
d <- dim(x)[2]
n <- dim(x)[1]
if(! method %in% c("Penalized","HJPS2","HJPSorig")) method <- "Penalized"
psiofx <- switch(basis,Quadratic=svd(cbind(rep(1,n),x,t(apply(x,1,expand))))$u,
                         Linear=svd(cbind(rep(1,n),x))$u,
			 svd(cbind(rep(1,n),x,t(apply(x,1,expand))))$u)
if(is.null(h0)) #h0<-.85*(d/n*prod(apply(apply(x,2,range),2,diff)/2))^(1/d)*sqrt(d)
# bandwidth in last line was approriate for uniform design only
# corrected to also give a better guess for "Gaussian" design 
#h0<-.85*(d/n*prod(sqrt(3*apply(x,2,var))))^(1/d)*sqrt(d)
# may be we need also something that works in case of heavy tailes
h0<-.85*(d/n*prod(apply(x,2,IQR)))^(1/d)*sqrt(d)
if(show>m) show<-m
if(graph && !is.null(fx)) oldpar <- par(mfrow=c(1,2)) 
if(is.null(cw)) cw <- max(6,log(n))
if(cw>=1) cw <- 1/cw
if(length(y)!=n) stop("wrong number of observations")
w0 <- sqrt(c(0.1252,0.1107,0.1,0.09075,0.08332,0.07703,
          0.07119,0.06671,0.06273,0.05868,0.05559))
if(d>(length(w0)+1)) w0d <- sqrt(1/(d+6))
else w0d <- w0[d-1]
h <- h0
hmax <- max(h0*n^(1/max(4,d)),h0*ch^100)
#  ch^100 is about a_h^100, allows to use step 4 for an increase of h about 100 times
#  this is just to allow use a very small initial guess for h_0, i.e. step 4 guarantees
#  that h_0 is increased until enough points with positive weights are available before starting
#  iterations
rho <- rho0
rhomin <- rho0*(dim(x)[1])^(-1/3)*crhomin
crho <- .846
cumlam <- NULL
nmean <- NULL
if(is.null(psiofx)) psiofx <- svd(cbind(rep(1,n),x),nv=0)$u
L <- dim(psiofx)[2]
bhat <- matrix(0,d,L)
ypshat <- matrix(0,d,d)
if(method=="Penalized"){
Mhat <- diag(L)
Pik <-cm*rho^2*diag(L)+Mhat
}
nw0<-n*w0d*cw
while((rho>rhomin&h<hmax)){
#  h>hmax  should occur much later than rho<rhomin, it is tested here only to 
#  stop the algorithm if the strategy using step 4 fails.
z <- edrk(x,y,n,d,psiofx,rho,h,crho,ch,bhat,ypshat,nw0,method,trace)
rho <- z$rho
h <- z$h
bhat <- z$bhat
nw0 <- z$nw0
nmean <- c(nmean,z$nmean)
lam  <- svd(bhat)$d
if(method=="Penalized"){
bhatPik <- bhat%*%Pik
svdbhatPik <- svd(bhatPik)
ypshat <- svdbhatPik$u%*%diag(svdbhatPik$d^2)%*%t(svdbhatPik$u)/svdbhatPik$d[1]^2
if(m>1){
Mhat <- svdbhatPik$v[,1:m]%*%t(svdbhatPik$v[,1:m])
} else {
Mhat <- outer(svdbhatPik$v[,1],svdbhatPik$v[,1],"*")
}
Pik <- cm*rho^2*diag(L)+Mhat
bhat <- bhat/max(lam)
} else {
bhat <- bhat/max(lam)
ypshat <- bhat%*%t(bhat)
}
lam  <- lam/max(lam)
if(graph) {
zg<-list(x=x,y=y,fhat=z$fhat,bhat=bhat)
    plot.edr(zg,show,ylab="Y",title=paste("fhat(x); Bandwidth= ",signif(h/ch,3)," Rho=",signif(rho/crho,3)))
    if(!is.null(fx)){
       zg$y <- fx
       zg$fhat <- fx
    plot.edr(zg,show,ylab="True f(x)",title=paste("True f(x): Bandwidth= ",signif(h/ch,3)," Rho=",signif(rho/crho,3)),sm=FALSE)
       }
    }
if(trace) {
cat("h=",signif(h/ch,3)," rho=",signif(rho/crho,3)," cummulative normalized EV",signif(cumsum(lam)[1:m]/sum(lam),3))
if(!is.null(R)) {
cat(" Risk 1:",signif(loss.edr(edr.R(bhat,show),R,dim(R)[2]),3)," Risk 2:",signif(loss2.edr(edr.R(bhat,m),R,show,x),3))
}
cat("\n")
gc()
}
}
if(graph && !is.null(fx)) par(oldpar) 
z<-list(x=x,y=y,bhat=bhat,fhat=z$fhat,cumlam=cumsum(lam)[1:m]/sum(lam),nmean=nmean,h=h/ch,rho=rho/crho,
h0=h0,rho0=rho0,cm=cm,call=args)
class(z)<-"edr"
z
}


edrk <- function(x,y,n,d,psiofx,rho,h,crho,ch,bhat,ypshat,nw0,method,trace){
#
#  one iteration of the algorithm
#
   eps <- .0000001
   d1 <- d+1
   L <- dim(psiofx)[2]
   sk2 <- diag(d)+ypshat/rho^2
   sk <- svd(sk2,nv=0)
   sk <- sk$u%*%diag(sqrt(sk$d))%*%t(sk$u)
#
#  step 3 and 4
#
   wrongh <- TRUE
   counthincrease <- 0
   wij <- x%*%sk
   Kksi <- dist(wij)^2
   fhat <- y
   while(wrongh){
   if(method=="HJPSorig"){
   wi<-.Fortran("edrstp3a",as.integer(d),
                           as.integer(n),
			   as.integer(d1),
                           as.double(wij),
			   as.double(Kksi/h^2),
			   double(n),
		           double(d1*n),
			   wi=double(n),
			   double(d1),
			   double(10*(3*d1+max(n,6*d1))),
		           integer(8*d1),
			   PACKAGE="EDR")$wi
   } else {
   wi<-.Fortran("edrstp3a",as.integer(d),
                           as.integer(n),
			   as.integer(d1),
                           as.double(x),
			   as.double(Kksi/h^2),
			   double(n),
		           double(d1*n),
			   wi=double(n),
			   double(d1),
			   double(10*(3*d1+max(n,6*d1))),
		           integer(8*d1),
			   PACKAGE="EDR")$wi
    }
    swi <- sum(wi)
    if(swi>=nw0) wrongh <- FALSE
    else {
       h <- h*ch
       if(trace) cat("increase bandwidth to",h,"swi,nw0",swi,nw0,"\n")
       counthincrease <- counthincrease+1
       if(counthincrease>(d+2)) wrongh <- FALSE  
       }
    }
   fx <- matrix(0,n,d)
   lll <- 0
   Kksi <- Kksi/h^2
   z<-.Fortran("edrstp3b",as.integer(d),
                          as.integer(n),
			  as.integer(d1),
                          as.double(x),#wij
			  as.double(Kksi),#kksi
			  as.double(y),#y
			  double(n),#kksii
		          double(d1*n),#mat
			  double(d1),#s
			  double(d1*d1),#u
			  double(d1*n),#v
			  double(5*(3*d1*d1+max(n,4*d1*d1+4*d1))),#work
		          integer(8*d1),#iwork
			  fx=as.double(fx),
			  fhat=as.double(fhat), 
			  lll=as.double(lll),
			  yw=double(n),#yw
			  PACKAGE="EDR")[c("fx","fhat","lll")]
      fx<-matrix(z$fx,n,d)
      lll<-z$lll
   for(i in 1:L){
      bhat[,i] <- rep(1,n)%*%(fx*psiofx[,i]*wi)/swi
#      apply(fx1*psiofx[,i]*wi,2,sum)/swi
      }
   z<-list(rho=rho*crho,h=h*ch,bhat=bhat,fhat=z$fhat,nmean=lll/n,nw0=nw0*crho)
   }

plot.edr <- function(x,m=1,ylab="Y",title="",sm=require(sm),...){
    xx<-x$x
    yy<-x$y
    n <- length(yy)
    bhat<-x$bhat
    fhat <- x$fhat
    if(!sm){
    if(m>=2 & !require(akima, quietly = TRUE)){
        warning(" Package sm or akima needed for 2D interpolation \n
             Note that Package akima is distributed under ACM license \n(see
             http://www.acm.org/publications/policies/softwarecrnotice.is )\n
             Using univariate projections instead of 2D\n")
        m <- 1
    }
    if(m==1){
    plot(xx%*%t(edr.R(bhat,1)),yy,xlab="R_m%*%x",ylab=ylab,...)
    xxx<-xx%*%t(edr.R(bhat,1))
    points(xxx,fhat,col=2)
    }
    if(m>=2){
     xr<-xx%*%t(edr.R(bhat,2))
     ngrid<-max(100,length(yy)^(2/3))
     xo<-seq(min(xr[,1]),max(xr[,1]),length=ngrid)
     yo<-seq(min(xr[,2]),max(xr[,2]),length=ngrid)
     z<-interp(xr[,1],xr[,2],fhat,xo=xo,yo=yo)
     image(z,xlab="First projection",ylab="Second projection",...)
     contour(z,add=TRUE)
    }
    } else {
    if(m==1) ngrid<-length(yy) else ngrid<-max(100,n^(2/3))
    z <- sm.regression(xx%*%t(edr.R(bhat,m)),yy,h=h.select(xx%*%t(edr.R(bhat,m)),yy,method="cv"),display="none",ngrid=ngrid)
    if(m==1){
    plot(xx%*%t(edr.R(bhat,m)),yy,xlab="R_m%*%x",ylab=ylab,...)
    lines(z$eval.points,z$estimate,col=2)
    }
    if(m>=2){
     image(z$estimate,xlab="First projection",ylab="Second projection",...)
     contour(z$estimate,add=TRUE)
    }
    }
    title(title)
invisible(NULL)
}

summary.edr <- function(object,m=1,R=NULL,...){
cat("Basis of the Effective Dimension Reduction Space, m=",m,"\n\n")
Rhat <- edr.R(object,m)
print(signif(Rhat,4))
if(!is.null(R)){
cat("Distance ||R %*% (I- hat{P}_m)||/||R|| between true and estimated EDR (m=",m,") \n\n")
loss1 <- loss.edr(edr.R(object,m),R,dim(R)[2])
cat(signif(loss1,3),"\n\n")
cat("Distance defined by Li (1992) between true and estimated EDR (m=",m,") \n\n")
loss2 <- loss2.edr(edr.R(object,m),R,m,object$x)
cat(signif(loss2,3),"\n\n")
} else {
loss1 <- loss2 <- NULL
}
cat("Cumulative sum of normalized eigenvalues of Bhat:",object$cumlam,"\n\n")
if(!is.null(object$cvmse)){
cat("Estimated MSEP",signif(object$cvmse,3),"\n\n")
cat("Estimated MAEP",signif(object$cvmae,3),"\n\n")
if(!is.null(object$hsmopt)) cat("Bandwidth used within EDR",signif(object$hsmopt,3),"\n\n")
}
list(Rhat=Rhat,cumlam=object$cumlam,loss1=loss1,loss2=loss2)
}

print.edr <- function(x,m=1,R=NULL,...){
cat("Basis of the Effective Dimension Reduction Space, m=",m,"\n\n")
print(signif(edr.R(x,m),4))
if(!is.null(R)){
cat("Distance ||R %*% (I- hat{P}_m)||/||R|| between true and estimated EDR (m=",m,") \n\n")
cat(signif(loss.edr(edr.R(x,m),R,dim(R)[2]),3),"\n\n")
cat("Distance defined by Li (1992) between true and estimated EDR (m=",m,") \n\n")
cat(signif(loss2.edr(edr.R(x,m),R,m,x$x),3),"\n\n")
}
cat("Cumulative sum of normalized eigenvalues of Bhat:",signif(x$cumlam,3),"\n\n")
cat("Initial bandwidth:",signif(x$h0,3),"\n\n")
cat("Final bandwidth:",signif(x$h,3),"\n\n")
cat("Initial value of rho:",signif(x$rho0,3),"\n\n")
cat("Final value of rho:",signif(x$rho,3),"\n\n")
if(!is.null(x$cvmse)){
cat("Estimated MSEP",signif(x$cvmse,3),"\n\n")
cat("Estimated MAEP",signif(x$cvmae,3),"\n\n")
if(!is.null(x$hsmopt)) cat("Bandwidth used within EDR",signif(x$hsmopt,3),"\n\n")
}
invisible(NULL)
}

edr.R <- function(B,m){
# compute matrix R_m such that  R_m%*%x lies in the m-dimensional
# subspace 
if(class(B)=="edr") B <- B$bhat
#
#   extract matrix Bhat if B is an object returned by function edr 
#   otherwise assume that B is the matrix Bhat
#
dimbhat <- dim(B)
if(is.null(dimbhat)||length(dimbhat)!=2) return(warning("Argument B should be of class edr or a matrix"))
if(dimbhat[2]<m) return(warning("Argument m is to large"))
t(svd(B,nv=0)$u[,1:m,drop=FALSE])
}

loss.edr <- function(Rstar,R,d){
#  compute ||R %*% (I- hat{P}_m)||/||R|| with P_m being the projection matrix  
#  for the estimated EDR
#  this defines a distance between spaces spanned by Rstar and the estimated 
#  EDR
zz <- svd(t(Rstar),nv=0)$u
(sum(svd(R%*%(diag(d)-zz%*%t(zz)))$d^2)/sum(svd(R)$d^2))
}
loss2.edr <- function(Rstar,R,m,x){
# thats cancor following Li (1992) 
1-mean(cancor(x%*%t(Rstar[1:m,,drop=FALSE]),x%*%t(R[1:m,,drop=FALSE]))$cor^2)
}
