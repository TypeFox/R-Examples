chi.int2.p <-
function(p,a,c1)
return( -exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*dchisq(c1^2,p+a)*2*c1 )

