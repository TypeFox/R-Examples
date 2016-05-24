pgfwaring <-
function(s,params) {
    require(hypergeo)
k<-s[abs(s)>1]


if (length(k)>0)

warning("At least one element of the vector s are out of interval [-1,1]")

if (length(params)<2) stop("At least one value in params is missing")
if (length(params)>2) stop("The length of params is 2")

    cc<-params[1]
    a<-params[2]
if (cc<=0)
stop ("Parameter c must be positive")
if (a<=0)
stop ("Parameter a must be positive")


    (cc-a)/cc*Re(hypergeo(1,a,cc+1,s))
}
