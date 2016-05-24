calculate.mle<-function(d, model="poisson")
{
	mle<-switch(model,
	multinomial=list(n=sum(d$Observed), p=scale(d$Expected, center=FALSE, scale=sum(d$Expected))),
	poisson=list(n=length(d$Observed), lambda=d$Expected)
	)
	
	if(model=="negbin")
	{
		smth<-empbaysmooth(d$Observed, d$Expected)
		mle<-list(n=length(d$Observed), nu=smth$nu, alpha=smth$alpha,
		size=smth$nu,  prob=smth$alpha/(smth$alpha+d$Expected) ) 
	}

	return(mle)
}
