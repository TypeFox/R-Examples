pgfthomas <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<2) stop("At least one value in params is missing")
if (length(params)>2) stop("The length of params is 2")

    lambda<-params[1]
    theta<-params[2]
if (lambda<=0)
stop ("Parameter lambda must be positive")
if (theta<=0)
stop ("Parameter theta must be positive")

    exp(lambda*(s*exp(theta*(s-1))-1))
}
