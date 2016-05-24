require(R2Cuba)

# Definitions of Genz test functions

f1 <- function(x,phw=NULL,cpar,wpar)
  cos(sum(x*cpar)+2*pi*wpar[1])
f2 <- function(x,cpar,wpar,...)
  1/prod((x-wpar)^2+1/cpar^2)
f3 <- function(x,cpar,wpar,...)
  1/(1+sum(x*cpar))^(length(x)+1)
f4 <- function(x,cpar,wpar,...)
  exp(-sum(cpar^2*(x-wpar)^2))
f5 <- function(x,cpar,wpar,...)
  exp(-sum(cpar*abs(x-wpar)))
f6 <- function(x,cpar,wpar,...)
  exp(sum(x*cpar))*prod(x<wpar)

# All Genz functions have known expressions for their
# integrals

allA <- function(d) { # utility function
  # compute all subsets from {1,...,d}
  ival <- vector("list",length=d)
  for(i in 1:d) ival[[i]] <- 0:1
  names(ival) <- paste("i",1:d,sep="")
  as.matrix(expand.grid(ival))
}

int1 <- function(d,cpar,wpar) {
  int1Contrib <- function(A,cpar,wpar) {
    cosines <- ifelse(cpar!=0,(1-cos(cpar))/cpar,0)
    sines <- ifelse(cpar!=0,sin(cpar)/cpar,1)
    cos((sum(A)%%4)*pi/2+2*pi*wpar[1])*prod(ifelse(A,cosines,sines))
  }
  sum(apply(allA(d),1,int1Contrib,cpar,wpar[1]))
}

int2 <- function(d,cpar,wpar) {
  prod(cpar)*prod(atan(cpar*(1-wpar))+atan(cpar*wpar))
}

int3 <- function(d,cpar,wpar) {
  buf <- (-1)^d/(prod(1:d)*prod(cpar))
  contribA <- function(A,cpar,wpar) {
    prod(rep(-1,sum(A)))/(1+sum(cpar[!A]))
  }
  buf*sum(apply(allA(d),1,contribA,cpar=cpar,wpar=wpar))
}

int4 <- function(d,cpar,wpar) {
  buf <- pi^(d/2)/prod(cpar)
  buf*prod(pnorm(sqrt(2)*cpar*(1-wpar))-pnorm(-sqrt(2)*cpar*wpar))
}

int5 <- function(d,cpar,wpar) {
 u <- pmin(pmax(0,wpar),1)
 prod((2*cosh(cpar*(u-wpar))-exp(-cpar*wpar)-exp(cpar*(wpar-1)))/cpar)
}

int6 <- function(d,cpar,wpar) {
  u <- pmin(pmax(0,wpar),1)
  prod((exp(cpar*u)-1)/cpar)
}

## choose some random parameters for the test functions

cpar <- runif(2)
wpar <- runif(2)
cpar <- cpar/sum(cpar)
cpar <- 2*cpar

## first test function (see plot)

xgrid <- as.matrix(expand.grid(x1=seq(0,1,length=100),x2= seq(0,1,length=100)))
fval <- apply(xgrid,1,f1,cpar=cpar,wpar=wpar)
fval <- matrix(fval,100,100)
par(mfrow=c(2,3))
image(z=fval,asp=1,main=paste("f1 integral =",signif(int1(2,cpar,wpar),4)))

## Numerical integration

cuhre(  ndim=2,ncomp=1,integrand=f1,cpar=cpar,wpar=wpar,flags=list(verbose=0))
vegas(  ndim=2,ncomp=1,integrand=f1,cpar=cpar,wpar=wpar,flags=list(verbose=0))
suave(  ndim=2,ncomp=1,integrand=f1,cpar=cpar,wpar=wpar,flags=list(verbose=0))
divonne(ndim=2,ncomp=1,integrand=f1,cpar=cpar,wpar=wpar,flags=list(verbose=0))

## Second test function

fval <- apply(xgrid,1,f2,cpar=cpar,wpar=wpar)
fval <- matrix(fval,100,100)
image(z=fval,asp=1,main=paste("f2 integral =",signif(int2(2,cpar,wpar),4)))

## Numerical integration

cuhre(  ndim=2,ncomp=1,integrand=f2,cpar=cpar,wpar=wpar,flags=list(verbose=0))
vegas(  ndim=2,ncomp=1,integrand=f2,cpar=cpar,wpar=wpar,flags=list(verbose=0))
suave(  ndim=2,ncomp=1,integrand=f2,cpar=cpar,wpar=wpar,flags=list(verbose=0))
divonne(ndim=2,ncomp=1,integrand=f2,cpar=cpar,wpar=wpar,flags=list(verbose=0))

## Third test function

fval <- apply(xgrid,1,f3,cpar=cpar,wpar=wpar)
fval <- matrix(fval,100,100)
image(z=fval,asp=1,main=paste("f3 integral =",signif(int3(2,cpar,wpar),4)))

## Numerical integration

cuhre(  ndim=2,ncomp=1,integrand=f3,cpar=cpar,wpar=wpar,flags=list(verbose=0))
vegas(  ndim=2,ncomp=1,integrand=f3,cpar=cpar,wpar=wpar,flags=list(verbose=0))
suave(  ndim=2,ncomp=1,integrand=f3,cpar=cpar,wpar=wpar,flags=list(verbose=0))
divonne(ndim=2,ncomp=1,integrand=f3,cpar=cpar,wpar=wpar,flags=list(verbose=0))

## Fourth test function

fval <- apply(xgrid,1,f4,cpar=cpar,wpar=wpar)
fval <- matrix(fval,100,100)
image(z=fval,asp=1,main=paste("f4 integral =",signif(int4(2,cpar,wpar),4)))

## Numerical integration

cuhre(  ndim=2,ncomp=1,integrand=f4,cpar=cpar,wpar=wpar,flags=list(verbose=0))
vegas(  ndim=2,ncomp=1,integrand=f4,cpar=cpar,wpar=wpar,flags=list(verbose=0))
suave(  ndim=2,ncomp=1,integrand=f4,cpar=cpar,wpar=wpar,flags=list(verbose=0))
divonne(ndim=2,ncomp=1,integrand=f4,cpar=cpar,wpar=wpar,flags=list(verbose=0))

## Fifth test function

fval <- apply(xgrid,1,f5,cpar=cpar,wpar=wpar)
fval <- matrix(fval,100,100)
image(z=fval,asp=1,main=paste("f5 integral =",signif(int5(2,cpar,wpar),4)))

## Numerical integration

cuhre(  ndim=2,ncomp=1,integrand=f5,cpar=cpar,wpar=wpar,flags=list(verbose=0))
vegas(  ndim=2,ncomp=1,integrand=f5,cpar=cpar,wpar=wpar,flags=list(verbose=0))
suave(  ndim=2,ncomp=1,integrand=f5,cpar=cpar,wpar=wpar,flags=list(verbose=0))
divonne(ndim=2,ncomp=1,integrand=f5,cpar=cpar,wpar=wpar,flags=list(verbose=0))

## Sixth test function

fval <- apply(xgrid,1,f6,cpar=cpar,wpar=wpar)
fval <- matrix(fval,100,100)
image(z=fval,asp=1,main=paste("f6 integral =",signif(int6(2,cpar,wpar),4)))

## Numerical integration

cuhre(  ndim=2,ncomp=1,integrand=f6,cpar=cpar,wpar=wpar,flags=list(verbose=0))
vegas(  ndim=2,ncomp=1,integrand=f6,cpar=cpar,wpar=wpar,flags=list(verbose=0))
suave(  ndim=2,ncomp=1,integrand=f6,cpar=cpar,wpar=wpar,flags=list(verbose=0))
divonne(ndim=2,ncomp=1,integrand=f6,cpar=cpar,wpar=wpar,flags=list(verbose=0))
