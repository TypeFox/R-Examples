negbin.sim<-function(data, mle=NULL)
{
	if(is.null(mle))
	{
		mle<-empbaysmooth(data$Observed, data$Expected)
		mle<-list(n=mle$n, size=mle$nu, 
			prob=mle$alpha/(mle$alpha+data$Expected) )
	}
	
	data$Observed<-rnbinom(n=mle$n, size=mle$size, prob=mle$prob)
	return(data)
}

