pgfhypergeometric <-
function(s,params) {
    require(hypergeo)

k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<3) stop("At least one value in params is missing")
if (length(params)>3) stop("The length of params is 3")


    m<-params[1]
    n<-params[2]
    p<-params[3]

if (m<0)
     stop("Parameter m must be positive")
 if(!(abs(m-round(m))<.Machine$double.eps^0.5))
stop("Parameter m must be positive integer")

if (n<0)
     stop("Parameter n must be positive")
 if(!(abs(n-round(n))<.Machine$double.eps^0.5))
stop("Parameter n must be positive integer")
if ((p>=1)|(p<=0))
stop ("Parameter p belongs to the interval (0,1)")

if (m<n)
stop ("Parameter m is greater or equal then n ")
    Re(hypergeo(-n,-m*p,-m,1-s))
}
