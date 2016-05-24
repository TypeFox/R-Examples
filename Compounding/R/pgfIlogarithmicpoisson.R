pgfIlogarithmicpoisson <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<2) stop("At least one value in params is missing")
if (length(params)>2) stop("The length of params is 2")




    theta<-params[1]
    lambda<-params[2]
if ((theta>=1)|(theta<=0))
stop ("Parameter theta belongs to the interval (0,1)")

if (lambda<=0)
stop ("Parameter lambda must be positive")

    zval<-(1-theta^s)/(1-theta)
    1+log(zval)/lambda
}
