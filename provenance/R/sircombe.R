# Sircombe and Hazelton:
# Density estimator with sample-point adaptive bandwidth
adapt.f <- function(x,h,eval){
    f <- eval*0
    for (i in 1:length(eval)){f[i] <- mean(stats::dnorm(x-eval[i],0,sd=h))}
    return(f)
}

# Sircombe and Hazelton:
sig2.con <- function(x,sigx,UCV=TRUE){
    sig2max <- max(sigx^2)
    h <- 1.06*min(stats::sd(x),stats::IQR(x)/1.34)/length(x)^0.2
    if (UCV) h <- stats::bw.ucv(x,lower=h/20,upper=h)
    sig2.con <- sig2max+h^2
    return(sig2.con)
}

# Sircombe and Hazelton:
# Calculate integrated squared KDEs
Rf <- function(x,sig2,c.con){
    h1 <- sqrt(c(outer(c.con-sig2,c.con-sig2,"+")))
    xdiff <- c(outer(x,x,"-"))
    rf <- mean(stats::dnorm(xdiff,sd=h1))
    return(rf)
}

# Sircombe and Hazelton:
# calculate 'c' value for data (this returns c^2)
getc2 <- function(x){
    n <- length(x$x)
    c2 <- 0
    for (i in 1:n){
        foo <- sig2.con(x$x[[i]],x$err[[i]])
        c2 <- max(c2,foo)
    }
    return(c2)
}

#' Sircombe and Hazelton distance
#'
#' Calculates Sircombe and Hazelton's L2 distance between the Kernel
#' Functional Estimates (KFEs, not to be confused with Kernel Density
#' Estimates!) of two samples with specified analytical uncertainties
#'
#' @param x an object of class \code{distributional}
#' @param i index of the first sample
#' @param j index of the second sample
#' @param c.con smoothing bandwidth of the kernel functional estimate
#' @return a scalar value expressing the L2 distance between the KFEs
#' of samples i and j
#' @author Keith Sircombe and Martin Hazelton
#' @references Sircombe, K. N., and M. L. Hazelton. "Comparison of
#' detrital zircon age distributions by kernel functional estimation."
#' Sedimentary Geology 171.1 (2004): 91-111.
#' @examples
#' datfile <- system.file("DZ.csv",package="provenance")
#' errfile <- system.file("DZerr.csv",package="provenance")
#' DZ <- read.distributional(datfile,errfile)
#' d <- SH.diss(DZ,1,2)
#' print(d)
#' @seealso KS.diss
#' @export
SH.diss <- function(x,i,j,c.con=0){
    X <- x$x[[i]]
    Y <- x$x[[j]]
    sigX <- x$err[[i]]
    sigY <- x$err[[j]]
    sig2X <- sigX^2
    sig2Y <- sigY^2
    if (c.con <= 0) c.con <- max(sig2.con(X,sigX),sig2.con(Y,sigY))
    rfX <- Rf(X,sig2X,c.con)
    rfY <- Rf(Y,sig2Y,c.con)
    h1 <- sqrt(c(outer(c.con-sig2X,c.con-sig2Y,"+")))
    XYdiff <- c(outer(X,Y,"-"))
    Itmp <- mean(stats::dnorm(XYdiff,sd=h1))
    dXY <- rfX+rfY-2*Itmp
    dXY <- sqrt(dXY)
    return(dXY)
}
