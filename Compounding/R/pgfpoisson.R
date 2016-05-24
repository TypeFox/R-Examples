pgfpoisson <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)
warning("At least one element of the vector s are out of interval [-1,1]")
if (missing(params)) stop("Distribution parameters are not defined")
    theta<-params[1]
if (theta<=0)
stop ("Parameter of Poisson distribution must be positive")
    exp(theta*(s-1))
}
