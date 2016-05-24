rmix.norm <-
function (n,alpha,mu,sigma=rep(1,length(alpha))) 
#n:      sample size.
#alpha:  vector of mixture probabilities.
#mu:     vector of means of each component.
#sigma:  vector of standard deviation of each component.
{
	m=length(alpha)
	alpha=alpha/sum(alpha)
	data=c()

	nindex=rmultinom(1,n,alpha)

	for( i in 1:m)
		data=c(data,rnorm(nindex[i],mu[i],sigma[i]))
	data
}
