pmixexp <-
function(x,alpha,theta,alp0)
{
#This function computes the cdf of Exponential mixture minus a given probability alp0.
#
#x:   data, a vector of observed values.
#alpha:  vector of mixture probabilities.
#theta:  vector of means of each component.
#alp0:   a given probability.	
	sum(alpha*pexp(x,rate=1/theta))-alp0
}
