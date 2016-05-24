rYP <- function(th1, th2, n, aX) {
#### Please refer to the function simuDataYP() on the method
#### used here. We return  H^(-1) (Z), with Z been exp(1) r.v.
#### th1=exp(short term hazard ratio); th2=exp(long term hazard ratio)
#### th1 and th2 can either be of length n or 1. Should be th1=exp(beta1*Z),
####  th2=exp(beta2*Z), same Z. 


if( any(th1 <= 0) ) stop("th1 must all > 0")
if( any(th2 <= 0) ) stop("th2 must all > 0")
EaX2 <- exp( aX )
if(length(aX) != n) stop("check length of aX")
XUnif <- runif(n)
temp <- 1+(th2/th1)*( XUnif^(-1/(th2*EaX2)) -1)
## may also use the following two lines
## Xexp <- rexp(n)
## temp <- 1+( exp(Xexp/(th2*aX)) -1)* th2/th1
if (any(temp <= 0)) stop( "check inputs" )
return(log(temp))
}
