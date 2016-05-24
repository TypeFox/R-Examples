poisson.sim<-function(data, mle=NULL)
{
	if(is.null(mle))
		data$Observed<-rpois(length(data$Observed), lambda=data$Expected)
	else
		data$Observed<-rpois(n=mle$n, lambda=mle$lambda)
	return(data)
}
