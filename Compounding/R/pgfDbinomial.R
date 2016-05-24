pgfDbinomial <-
function(s,params) {

k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")
if (length(params)<2) stop("At least one value in params is missing")
if (length(params)>2) stop("The length of params is 2")

       n<-params[1]
       theta<-params[2]

if ((theta>=1)|(theta<=0))
stop ("Parameter theta belongs to the interval (0,1)")
if (n<0)
     stop("Parameter n must be positive")
 if(!(abs(n-round(n))<.Machine$double.eps^0.5))
stop("Parameter n must be positive integer")
    n*theta*(1-theta+theta*s)^(n-1)
}
