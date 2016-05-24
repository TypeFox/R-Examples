"pdiscp" <-
function(p,probs,n,s)
{
#
# PDISCP Predictive distribution of number of successes in future binomial
#  experiment with a discrete prior. PRED = PDISCP(P,PROBS,N,S) returns 
#  vector PRED of predictive probabilities, where P is the vector of 
#  values of the proportion, PROBS is the corresponding vector of
#  probabilities, N is the future binomial sample size, and S is the vector of
#  numbers of successes for which predictive probabilities will be computed.
#------------------------
# Written by Jim Albert
# albert@bgnet.bgsu.edu
# November 2004
#------------------------

pred=0*s;

for (i in 1:length(p))
{
 pred=pred+probs[i]*dbinom(s,n,p[i]);
}
return(pred)
}

