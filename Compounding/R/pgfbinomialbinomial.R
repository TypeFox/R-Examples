pgfbinomialbinomial <-
function(s,params) {
k<-s[abs(s)>1]


if (length(k)>0)
warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<4) stop("At least one value in params is missing")
if (length(params)>4) stop("The length of params is 4")


    p1<-params[1]
    p2<-params[2]
    m<-params[3]
    n<-params[4]
if ((p1>=1)|(p1<=0))
stop ("Parameter p1 belongs to the interval (0,1)")
if ((p2>=1)|(p2<=0))
stop ("Parameter p2 belongs to the interval (0,1)")
if (m<0)
     stop("Parameter m must be positive")
 if(!(abs(m-round(m))<.Machine$double.eps^0.5))
stop("Parameter m must be positive")
if (n<0)
     stop("Parameter n must be positive")
 if(!(abs(n-round(n))<.Machine$double.eps^0.5))
stop("Parameter n must be positive integer")
    (1-p1+p1*(1-p2+p2*s)^n)^m
}
