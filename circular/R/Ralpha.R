Ralpha <- function(x, n, alpha) {
  #x is the C part of the observed R (resultant length)
   if (n==1) stop('We are not able to provide sensible results for n=1')
   I2n <- function(x, n, lower=NULL, upper=NULL) {
     if (is.null(lower))
       lower <- 0
     if (is.null(upper)) {
       if (n < 1500)
         upper <-  100+10000/n
       else
         upper <- 50+10000/n
     }
     #x is R here while in the next x is the variable of integration
     temp <- function(x, n, lower, upper) {
       f1 <- integrate(function(x, R, n) besselJ(x, 0)^n*besselJ(R*x, 0)*x, lower=lower, upper=upper-2*0.99*log(n), R=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f2 <- integrate(function(x, R, n) besselJ(x, 0)^n*besselJ(R*x, 0)*x, lower=lower, upper=upper-0.99*log(n), R=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f3 <- integrate(function(x, R, n) besselJ(x, 0)^n*besselJ(R*x, 0)*x, lower=lower, upper=upper, R=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f4 <- integrate(function(x, R, n) besselJ(x, 0)^n*besselJ(R*x, 0)*x, lower=lower, upper=upper+0.99*log(n), R=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f5 <- integrate(function(x, R, n) besselJ(x, 0)^n*besselJ(R*x, 0)*x, lower=lower, upper=upper+2*0.99*log(n), R=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       median(c(f1,f2,f3,f4,f5))
     }
     sapply(X=x, FUN=temp, n=n, lower=lower, upper=upper)
   }
   I3n <- function(x, n, lower=NULL, upper=NULL) {
     if (is.null(lower))
       lower <- 0
     if (is.null(upper)) {
####       upper <- approx(x=c(4, 5, 10, 13, 15, 20, 28, 50, 60, 75, 100, 250, 500, 1000, 2000, 5000), y=c(6000, 2000, 1250, 1000, 700, 600, 500, 400, 350, 300, 250, 200, 150, 90, 70, 50), xout=n, method='constant', yleft=6000, yright=40, rule=2, f=1)$y
####       lma <- lm(I(y[-(1:2)]~I(1/x[-(1:2)]))
####       upper <-  114.3+10856.2/n
       if (n < 1500)
         upper <-  100+10000/n
       else
         upper <- 50+10000/n
     }
     #x is the C part here while in the next x is the variable of integration
     temp <- function(x, n, lower, upper) {
       f1 <- integrate(function(x, C, n) besselJ(x, 0)^n*cos(C*x), lower=lower, upper=upper-2*0.99*log(n+1), C=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f2 <- integrate(function(x, C, n) besselJ(x, 0)^n*cos(C*x), lower=lower, upper=upper-0.99*log(n+1), C=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f3 <- integrate(function(x, C, n) besselJ(x, 0)^n*cos(C*x), lower=lower, upper=upper, C=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f4 <- integrate(function(x, C, n) besselJ(x, 0)^n*cos(C*x), lower=lower, upper=upper+0.99*log(n+1), C=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       f5 <- integrate(function(x, C, n) besselJ(x, 0)^n*cos(C*x), lower=lower, upper=upper+2*0.99*log(n+1), C=x, n=n, subdivisions=4000, stop.on.error=FALSE)$value
       median(c(f1,f2,f3,f4,f5))
     }
     sapply(X=x, temp, n=n, lower=lower, upper=upper)
   }   
   lhs <- function(x, C, n) { #left hand side of the equation
     #x is the R0 while in the next x is the variable of integration
     temp <- function(x, C, n) integrate(function(x, C, n) x/sqrt(x^2-C^2)*I2n(x, n), lower=x, upper=n, C=C, n=n, subdivisions=4000, stop.on.error=FALSE)$value
     sapply(X=x, FUN=temp, C=C, n=n) 
   }
   equat <- function(x, C, n, alphaI3n) { #equation
      lhs(x=x, C=C, n=n) - alphaI3n
   }
   temp <- function(x, n, alpha) {
     alphaI3n <- alpha*I3n(x=x, n=n)
     uniroot(equat, lower=x, upper=n, C=x, n=n, alphaI3n=alphaI3n)$root
   }
   sapply(X=x, FUN=temp, n=n, alpha=alpha)
}


Ralphaapprox <- function(x, n, alpha) {
  if (n<3) stop('We are not able to provide sensible results for n<3')
  temp <- function(x, n, alpha) {
    if (n >=15 & x > 0 & x < n/3) {
      y <- sqrt(x^2+qchisq(alpha, df=1, lower.tail=FALSE)*0.5*n) #3.2
    } else if (x > n/2 & x < 3*n/4) {
      ff <- qf(alpha, df1=2, df2=2*n-2, lower.tail=FALSE)
      y <- (ff*n+(n-1)*x)/(n+ff-1) #3.3
    } else if (x > 5/6*n) {
      ff <- qf(alpha, df1=1, df2=n-1, lower.tail=FALSE)
      y <- (ff*n+(n-1)*x)/(n+ff-1) #3.4
    } else {
      y <- NA
    }
    return(y)
  }
  sapply(X=x, FUN=temp, n=n, alpha=alpha)
}

