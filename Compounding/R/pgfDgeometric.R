pgfDgeometric <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)
warning("At least one element of the vector s are out of interval [-1,1]")


if (length(params)>1) stop("The length of params is 1")

    theta<-params[1]
if ((theta>=1)|(theta<=0))
stop ("Parameter of geometric distribution must belong  to the interval (0,1)")
    theta*(1-theta)/(1-(1-theta)*s)^2
}
