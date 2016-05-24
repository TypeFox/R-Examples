laplace <- function(dist, g, k){
if(dist=="gamma")  (1+g/k)^(-k)
else if(dist=="positives") exp(-g^k)
else if(dist=="power") 1-(1-exp(-g))^k
else if(dist=="lognormal") gh(g,0, k)/sqrt(2*pi)
else stop("Unrecognized frailty distribution")
}