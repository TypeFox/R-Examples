rmix.exp <-
function (n,alpha,theta) 
#n:      sample size.
#alpha:  vector of mixture probabilities.
#theta:  vector of parameters of each component.
{
	m=length(alpha)
	alpha=alpha/sum(alpha)
	data=c()

	nindex=rmultinom(1,n,alpha)

	for( i in 1:m)
		data=c(data,rexp(nindex[i],1/theta[i]))
	data
}
