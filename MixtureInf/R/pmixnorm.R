pmixnorm <-
function(x,alpha,theta,alp0)
{
#This function computes the cdf of Normal mixture minus a given probability alp0.
#
#x:   data, a vector of observed values.
#alpha:  vector of mixture probabilities.
#theta:  vector of means of each component.
#alp0:   a given probability. 	
	sum(alpha*pnorm(x,theta,1))-alp0
}
