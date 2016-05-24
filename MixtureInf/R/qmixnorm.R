qmixnorm <-
function(alpha,theta,alp0)
{
#This function computes the alp0-quantile of Normal mixture.
#
#alpha:  vector of mixture probabilities.
#theta:  vector of means of each component.
#alp0:   a given probability.		
	uniroot(pmixnorm,c(qnorm(alp0,min(theta)),qnorm(alp0,max(theta))),alpha=alpha,theta=theta,alp0=alp0)$root 
}
