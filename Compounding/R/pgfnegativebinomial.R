pgfnegativebinomial <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)
warning("Some  elements of the vector s are out of interval [-1,1]")

if (length(params)<2) stop("At least one value in params is missing")
if (length(params)>2) stop("The length of params is 2")
    theta<-params[1]
    k<-params[2]
if ((theta>=1)|(theta<=0))
stop ("Parameter theta belongs to the interval (0,1)")
if (k<=0)
stop("Parameter k must be positive")
    (theta/(1-(1-theta)*s))^k
}
