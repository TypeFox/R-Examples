pbetap=function(ab,n,s)
{
#
# PBETAP Predictive distribution of number of successes in future binomial
#	experiment with a beta prior. PRED = PBETAP(AB,N,S) returns a vector 
#	PRED of predictive probabilities, where AB is the vector of beta
#	parameters, N is the future binomial sample size, and S is the vector of
#	numbers of successes for which predictive probabilities will be computed.
#------------------------
# Written by Jim Albert
# albert@bgnet.bgsu.edu
# November 2004
#------------------------

pred=0*s;
a=ab[1]; b=ab[2];

lcon=lgamma(n+1)-lgamma(s+1)-lgamma(n-s+1);

pred=exp(lcon+lbeta(s+a,n-s+b)-lbeta(a,b));

return(pred)
}
