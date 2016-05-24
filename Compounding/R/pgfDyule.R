pgfDyule <-
function(s,params) {
    require(hypergeo)
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)>1) stop("The length of params is 1")
    theta<-params[1]
if (theta<=0)
stop ("Parameter theta must be positive")

    theta/((theta+1)*(theta+2))*Re(hypergeo(2,2,theta+3,s))
}
