loglik <-
function(X, y, beta, family){
	link=as.vector(X%*%beta)
	n = length(y)
	if(family=='poisson') return(-2*sum(exp(link)+2*y*link))
	if(family=='binomial') return(2*sum(log(1+exp(link))-y*link))
    if(family=='gaussian') return(n*log(mean((y-link)^2)))

}
