pgfDpoissonbinomial <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")
if (length(params)<3) stop("At least one value in params is missing")
if (length(params)>3) stop("The length of params is 3")
    theta<-params[1]
    p<-params[2]
    n<-params[3]
  
if (theta<=0)
stop ("Parameter theta must be positive")
if ((p>=1)|(p<=0))
stop ("Parameter p belongs to the interval (0,1)")
if (n<0)
     stop("Parameter n must be positive")
 if(!(abs(n-round(n))<.Machine$double.eps^0.5))
stop("Parameter n must be positive integer")
    n*theta*p*(1-p+p*s)^(n-1)*exp(theta*((1-p+p*s)^n-1))
}
