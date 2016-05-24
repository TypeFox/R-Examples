qmixexp <-
function(alpha,theta,alp0)
{
#This function computes the alp0-quantile of Exponential mixture.
#
#alpha:  vector of mixture probabilities.
#theta:  vector of means of each component.
#alo0:   a given probability.
	uniroot(pmixexp,c(0,qexp(alp0,rate=1/max(theta))),alpha=alpha,theta=theta,alp0=alp0)$root 
}
