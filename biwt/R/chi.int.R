chi.int <-
function(p,a,c1)
return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*pchisq(c1^2,p+a) )

