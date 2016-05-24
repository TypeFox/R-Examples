pgfpolyaaeppli <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<2) stop("At least one value in params is missing")
if (length(params)>2) stop("The length of params is 2")
    
theta<-params[1]
    p<-params[2]
if (theta<=0)
stop ("Parameter theta must be positive")

if ((p>=1)|(p<=0))
stop ("Parameter p belongs to the interval (0,1)")

    exp(theta/p*((1-p)/(1-p*s)-1))
}
