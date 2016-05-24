### R code from vignette source 'neldermead_manual.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: neldermead_manual.Rnw:47-49
###################################################
options(width=75)
require(neldermead)


###################################################
### code chunk number 2: neldermead_manual.Rnw:476-480 (eval = FALSE)
###################################################
##   myalgorithm <- function( this ){
##     ...
##     return(this)
##   }


###################################################
### code chunk number 3: neldermead_manual.Rnw:502-506 (eval = FALSE)
###################################################
##   mystoppingrule <- function( this , simplex ){
##   ...
##   return(list(this=this,terminate=terminate,status=status))
##   }


###################################################
### code chunk number 4: neldermead_manual.Rnw:572-588
###################################################
  quadratic <- function(x=NULL,index=NULL,fmsfundata=NULL){
    return(list(f=x[1]^2 + x[2]^2,
                g=c(),
                c=c(), 
                gc=c(),
                index=index,
                this=list(costfargument=fmsfundata)))
  }

  x0 <- transpose( c(1.0,1.0) )
  nm <- neldermead()
  nm <- neldermead.set(nm,'numberofvariables',2)
  nm <- neldermead.set(nm,'function',quadratic)
  nm <- neldermead.set(nm,'x0',x0)
  nm <- neldermead.search(nm)
  summary(nm)


###################################################
### code chunk number 5: neldermead_manual.Rnw:607-652
###################################################
  rosenbrock <- function(x=NULL,index=NULL,fmsfundata=NULL){
    return(list(f=100*(x[2]-x[1]^2)^2+(1-x[1])^2,
                g=c(),
                c=c(),
                gc=c(),
                index=index,
                this=list(costfargument=fmsfundata)))
  }
  x0 <- transpose(c(-1.2,1.0))
  nm <- neldermead()
  nm <- neldermead.set(nm,'numberofvariables',2)
  nm <- neldermead.set(nm,'function',rosenbrock)
  nm <- neldermead.set(nm,'x0',x0)
  nm <- neldermead.set(nm,'maxiter',200)
  nm <- neldermead.set(nm,'maxfunevals',300)
  nm <- neldermead.set(nm,'tolfunrelative',10*.Machine$double.eps)
  nm <- neldermead.set(nm,'tolxrelative',10*.Machine$double.eps)
  nm <- neldermead.set(nm,'simplex0method','axes')
  nm <- neldermead.set(nm,'simplex0length',1.0)
  nm <- neldermead.set(nm,'method','variable')
  nm <- neldermead.set(nm,'verbose',FALSE)
  nm <- neldermead.set(nm,'storehistory',TRUE)
  nm <- neldermead.set(nm,'verbosetermination',FALSE)
  nm <- neldermead.search(nm)
  
  xmin <- ymin <- -2.0 
  xmax <- ymax <- 2.0 
  nx <- ny <- 100
  stepy <- stepx <- (xmax - xmin)/nx
  ydata <- xdata <- seq(xmin,xmax,stepx)
  zdata <- apply(expand.grid(xdata,ydata),1,
                 function(x) neldermead.function(nm,transpose(x)))
  zdata <- matrix(zdata,ncol=length(ydata))
  optimpath <- matrix(unlist((neldermead.get(nm,'historyxopt'))),
                      nrow=2)
  optimpath <- data.frame(x=optimpath[1,],y=optimpath[2,])

  contour(xdata,ydata,zdata,levels=c(1,10,100,500,1000,2000))
  par(new=TRUE,ann=TRUE)
  plot(c(x0[1],optimpath$x[158]), c(x0[2],optimpath$y[158]),
       col=c('red','green'),pch=16,xlab='x[1]',ylab='x[2]',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax))
  par(new=TRUE,ann=FALSE)  
  plot(optimpath$x,optimpath$y,col='blue',type='l',
       xlim=c(xmin,xmax),ylim=c(ymin,ymax))


###################################################
### code chunk number 6: neldermead_manual.Rnw:694-716
###################################################
  quadratic <- function(x=NULL,index=NULL,fmsfundata=NULL){
    return(list(f=x[1]^2 + x[2]^2,
                g=c(),
                c=c(), 
                gc=c(),
                index=index,
                this=list(costfargument=fmsfundata)))
  }
  set.seed(0)
  x0 <- transpose(c(1.2,1.9))
  nm <- neldermead()
  nm <- neldermead.set(nm,'numberofvariables',2)
  nm <- neldermead.set(nm,'function',quadratic)
  nm <- neldermead.set(nm,'x0',x0)
  nm <- neldermead.set(nm,'verbose',FALSE)
  nm <- neldermead.set(nm,'storehistory',TRUE)
  nm <- neldermead.set(nm,'verbosetermination',FALSE)
  nm <- neldermead.set(nm,'method','box')
  nm <- neldermead.set(nm,'boundsmin',c(1,1))
  nm <- neldermead.set(nm,'boundsmax',c(2,2))
  nm <- neldermead.search(nm)
  summary(nm)


###################################################
### code chunk number 7: neldermead_manual.Rnw:750-784
###################################################
  michalewicz <- function(x=NULL,index=NULL,fmsfundata=NULL){
    f <- c()
    c <- c()
    if (index == 2 | index ==6) 
      f <- (x[1]-10)^3+(x[2]-20)^3
    
    if (index == 5 | index ==6)
      c <- c((x[1]-5)^2+(x[2]-5)^2 -100, 
          82.81-((x[1]-6)^2+(x[2]-5)^2))
    varargout <- list(f=f,
        g=c(),
        c=c, 
        gc=c(),
        index=index,
        this=list(costfargument=fmsfundata))
    return(varargout)
  }
  set.seed(0)
  x0 <- transpose(c(15,4.99))
  nm <- neldermead()
  nm <- neldermead.set(nm,'numberofvariables',2)
  nm <- neldermead.set(nm,'nbineqconst',2)
  nm <- neldermead.set(nm,'function',michalewicz)
  nm <- neldermead.set(nm,'x0',x0)
  nm <- neldermead.set(nm,'maxiter',300)
  nm <- neldermead.set(nm,'maxfunevals',1000)
  nm <- neldermead.set(nm,'simplex0method','randbounds')
  nm <- neldermead.set(nm,'boxnbpoints',3)
  nm <- neldermead.set(nm,'storehistory',TRUE)
  nm <- neldermead.set(nm,'method','box')
  nm <- neldermead.set(nm,'boundsmin',c(13,0))
  nm <- neldermead.set(nm,'boundsmax',c(20,10))
  nm <- neldermead.search(nm)
  summary(nm)


###################################################
### code chunk number 8: neldermead_manual.Rnw:804-833
###################################################
  negLL <- function(x=NULL, index=NULL, fmsfundata=NULL){
    mn <- x[1]
    sdv <- x[2]
    out <- -sum(dnorm(fmsfundata$data, mean=mn, sd=sdv, log=TRUE))
  
    return(list(f = out, 
           index = index,
           this=list(costfargument=fmsfundata)))
  }

  set.seed(12345)
  fmsfundata <- structure(
    list(data=rnorm(500,mean=50,sd=2)),
    class='optimbase.functionargs')

  x0 <- transpose(c(45,3))
  nm <- neldermead()
  nm <- neldermead.set(nm,'numberofvariables',2)
  nm <- neldermead.set(nm,'function',negLL)
  nm <- neldermead.set(nm,'x0',x0)
  nm <- neldermead.set(nm,'costfargument',fmsfundata)
  nm <- neldermead.set(nm,'maxiter',500)
  nm <- neldermead.set(nm,'maxfunevals',1500)
  nm <- neldermead.set(nm,'method','box')
  nm <- neldermead.set(nm,'storehistory',TRUE)
  nm <- neldermead.set(nm,'boundsmin',c(-100, 0))
  nm <- neldermead.set(nm,'boundsmax',c(100, 100))
  nm <- neldermead.search(this=nm)
  summary(nm)


###################################################
### code chunk number 9: neldermead_manual.Rnw:845-854
###################################################
  rosenbrock <- function(x=NULL){
    f <- 100*(x[2]-x[1]^2)^2+(1-x[1])^2
  }
  x0 <- c(-1.2,1.0)
  npts <- 6
  xmin <- c(-2,-2)
  xmax <- c(2,2)
  grid <- fmin.gridsearch(fun=rosenbrock,x0=x0,xmin=xmin,xmax=xmax,npts=npts,alpha=alpha)
  grid


