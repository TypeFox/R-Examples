pgfkattitypeh1 <-
function(s,params) {
    require(hypergeo)
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<3) stop("At least one value in params is missing")
if (length(params)>3) stop("The length of params is 3")


    theta<-params[1]
    a<-params[2]
    b<-params[3]

if (theta<=0)
stop ("Parameter theta must be positive")
if (a<=0)
stop ("Parameter a must be positive")
if (b<=0)
stop ("Parameter b must be positive")

    genhypergeo(a,a+b,theta*(s-1))
}
