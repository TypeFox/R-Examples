rmix.pois <-
function(n,alpha,theta)
#n:      sample size.
#alpha:  vector of mixture probabilities.
#theta:  vector of probabilities of success of each component.
{
	m=length(theta)
	alpha=alpha/sum(alpha)
	data=c()
	nindex=rmultinom(1,n,alpha)
	for(i in 1:m)
		data=c(data,rpois(nindex[i],theta[i]))
	data
}
