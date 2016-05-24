pgfDpascalpoisson <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<3) stop("At least one value in params is missing")
if (length(params)>3) stop("The length of params is 3")

theta<-params[1]
    mu<-params[2]
    a<-params[3]
if (theta<=0)
stop ("Parameter theta must be positive")
if (mu<=0)
stop ("Parameter mu must be positive")
if (a<=0)
stop ("Parameter a must be positive")
    mu*exp(theta*(s-1))*(1+mu/(a*theta)-mu/(a*theta)*exp(theta*(s-1)))^(-a-1)
}
