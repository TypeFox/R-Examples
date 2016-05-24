pgfDneymantypeb <-
function(s,params) {
    require(hypergeo)
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<2) stop("At least one value in params is missing")
if (length(params)>2) stop("The length of params is 2")
    theta<-params[1]
    lambda<-params[2]
if (theta<=0)
stop ("Parameter theta must be positive")
if (lambda<=0)
stop ("Parameter lambda must be positive")
    theta*lambda/2*genhypergeo(2,3,theta*(s-1))*pgfneymantypeb(s,params)
}
