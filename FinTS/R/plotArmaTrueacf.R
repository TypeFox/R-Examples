plotArmaTrueacf <- function(object, lag.max=20, pacf=FALSE, plot=TRUE,
                        xlab="lag", ylab=c("ACF", "PACF")[1+pacf],
                        ylim=c(-1, 1)*max(ACF), type="h",
                        complex.eps=1000*.Machine[["double.neg.eps"]],
                            ...){
##
##  1.  Extract AR part
##
  {
    if(is.list(object)){
      AR <- object$ar
      MA <- object$ma
    }
    else {
      AR <- object
      MA <- numeric(0)
    }
  }
##
##  2.  Compute roots and test stationarity 
##
  {
    if(length(AR)>0){
#     roots <- (1/polyroot(c(1, -AR)))
#     solve(polynomial(...)) is more accurate than polyroot
      library(polynom) 
      p <- polynomial(c(1, -AR))
      roots <- (1/solve(p))
# test stationarity
      if(any(abs(roots)>=1)){
        warning("root(s) outside the unit circle:  Process is NOT stationary.")
        out <- list(roots=roots, acf=NULL, periodicity=per0)
        if(pacf)names(out)[2] <- 'pacf' 
        return(out) 
      }
    }
    else
      roots <- complex(0)
  }
##
## 2.  Generate ACF
##
  ACF <- ARMAacf(ar=AR, ma=MA, lag.max=lag.max, pacf=pacf)
  if(plot){
    indx <- {
      if(pacf)1:lag.max
      else 0:lag.max
    }
    plot(indx, ACF, xlab=xlab, ylab=ylab, ylim=ylim, type=type, ...)
    abline(h=0)
  }
##
## 3.  Test periodicity    
##
  conjRoots <- findConjugates(roots)
  {
    if(length(conjRoots)>0){
      d <- abs(conjRoots)
      per <- (2*pi/acos(Re(conjRoots)/d))
      per0 <- data.frame(damping=d, period=per)
    }
    else{
      per0 <- data.frame(damping=numeric(0), period=numeric(0))
    }
  }
##
## 4.  Done    
##
  out <- list(roots=roots, acf=ACF, periodicity=per0)
  if(pacf)names(out)[2] <- 'pacf' 
  return(out) 
}
