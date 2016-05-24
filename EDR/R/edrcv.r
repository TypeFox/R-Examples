edrcv <- function(x,y,m=2,rho0=1,h0=NULL,ch=exp(.5/max(4,(dim(x)[2]))),
    crhomin=1,cm=4,method="Penalized",fit="sm",basis="Quadratic",cw=NULL,graph=FALSE,
    show=1,trace=FALSE,seed=1,cvsize=1,m0=min(m,2),hsm=NULL){
#
#   first prepare random grouping for Cross-validation 
#
args <- match.call()      
n <- length(y)
ind <- groups(n,cvsize,seed)

#
#  Now the CV-loop
#
if(fit=="sm") require(sm)
if(is.null(hsm)) lh <- 1 else lh <- length(hsm)
cvres <- matrix(0,n,lh)
cat("Start of CV loop \n Progress:")
ngroups <- length(unique(ind))
for(i in 1:ngroups){
xtest <- x[ind==i,]
ytest <- y[ind==i]
edrestimate <- edr(x[ind!=i,],y[ind!=i],m=m,rho0=rho0,h0=h0,ch=ch,
    crhomin=crhomin,cm=cm,method=method,basis=basis,cw=cw,
    show=show,trace=trace)
for(j in 1:lh){
cvfits <- predict.edr(edrestimate,xtest,m=m0,hsm[j],method=fit)
cvres[ind==i,j] <- ytest - cvfits$fhat
}
cat(" ",signif(i/ngroups*100,2),"%")
}
cat("\n")
# determine the optimal hsm
mae <- apply(abs(cvres),2,mean)
mse <- apply(cvres^2,2,mean)
if(!is.null(hsm)) hsmopt <- hsm[mae==min(mae)] else hsmopt <- NULL
cat("Estimate EDR using the full dataset")
#
#  Now produce all estimates using the full dataset
#
edrestimate <- edr(x,y,m=m,rho0=rho0,h0=h0,ch=ch,
    crhomin=crhomin,cm=cm,method=method,basis=basis,cw=cw,graph=graph,
    show=show,trace=trace)
edrestimate$call <- args
edrestimate$fhat <- predict.edr(edrestimate,x,m=m0,hsmopt,method=fit)
edrestimate$cvres <- cvres
edrestimate$cvmseofh <- mse
edrestimate$cvmaeofh <- mae
edrestimate$cvmse <- if(lh>1) mse[hsm==hsmopt] else mse
edrestimate$cvmae <- if(lh>1) mae[hsm==hsmopt] else mae
edrestimate$hsm <- hsm
edrestimate$hsmopt <- hsmopt
class(edrestimate) <- "edr"
edrestimate
}

groups <- function(n,size,seed=1){
set.seed(seed)
ox <- sample(1:n,n)
ngroups <- n%/%size
ind <- numeric(n)
k<-1
while(k <= n){
for( i in 1:ngroups ){
ind[ox[k]] <- i
k <- k+1
}
}
ind
}

predict.edr <- function(object,xest,m=1,h=NULL,method="sm",...){
    if(method=="sm"&!require(sm)) return("Please install package sm")
    if(method=="sm"&m>2) {
       method<-"direct"
#       warning("sm only handles up to 2 covariates, method reset to direct")
       }
    xx<-object$x
    yy<-object$y
    n <- length(yy)
    nest <- dim(xest)[1]
    bhat<-object$bhat
    fhat <- object$fhat
    P <- t(edr.R(bhat,m))
    xx <- xx%*%P
    xest <- xest%*%P
    if(m==1) {
       dim(xx) <- c(n,1)
       dim(xest) <- c(nest,1)
    }
    if(method=="sm"){
    if(m==1) ngrid<-length(yy) else ngrid<-max(100,n^(2/3))
    z <- sm::sm.regression(xx,yy,h=if(is.null(h)) sm::h.select(xx,yy,method="cv") else rep(h,m),display="none",eval.points=xest,eval.grid=FALSE)
    z <- list(x=xest,fhat=z$estimate)
    } else {
#  direct estimation
      if(is.null(h)) h <- object$h*object$rho
      wij <- t(xx)
      wijest <- t(xest)
      Kksi <- matrix(0,n,nest)
      for(i in 1:nest) Kksi[,i] <- apply((wij-wijest[,i])^2,2,sum)
      Kksi <- Kksi/h^2
      m1 <- m+1
      lll <- 0
      fhat <-.Fortran("edrstp3c",as.integer(m),
                          as.integer(n),
                          as.integer(nest),
			  as.integer(m1),
                          as.double(xx),#wij
                          as.double(xest),#wij
			  as.double(Kksi),#kksi
			  as.double(yy),#y
		          double(m1*n),#mat
			  double(m1),#s
			  double(m1*m1),#u
			  double(m1*n),#v
			  double(5*(3*m1*m1+max(n,4*m1*m1+4*m1))),#work
		          integer(8*m1),#iwork
			  fhat=double(nest), 
			  yw=double(n),#yw
			  PACKAGE="EDR")$fhat
    z <- list(x=xest,fhat=fhat)
    }
z
}
