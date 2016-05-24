rmix.binom <-
function(n,alpha,theta,size)
#n:      sample size.
#alpha:  vector of mixture probabilities.
#theta:  vector of probabilities of success of each component.
#size:   number of trials.
{
	m=length(theta)
	alpha=alpha/sum(alpha)
	data=c()
	nindex=rmultinom(1,n,alpha)
	for(i in 1:m)
		data=c(data,rbinom(nindex[i],size,theta[i]))
	data
}
