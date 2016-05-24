IPE <-function(para, method=NULL)
{
## para is the fitting results of one experiment
## method is "Poisson" or "NB"

	#1. Check data, method and initialpara are consistent
	if (sum(method=="Poisson")+sum(method=="NB")==0)
	{
		stop('A method must be given, either "Poisson" or "NB"!', call.=FALSE)
	}

	## we need to use interpolation since the discrete distribution will give some larger mass on some particular value
	## we interpolate between all values below cut
	T=10000
	cut=100
	n=500
	sum1=0
	sum2=0
	if (method=="Poisson")
	K=4
	if (method=="NB")
	K=6
	for (i in para[K]:cut)
	{
		if (method=="Poisson")
		{
			X.bg<-approx(x=c(i, i+1),y=c(dpois(i,para[3]),dpois(i+1,para[3])), n=n)
			Step<-X.bg$x[2]- X.bg$x[1]
			sum2=sum2+sum(ppois(i-para[K], para[2])*X.bg$y)*Step 
		}
		if (method=="NB")
		{
			X.bg<-approx(x=c(i, i+1),y=c(dnbinom(i,para[5], ,para[4]),dnbinom(i+1,para[5], ,para[4])), n=n)
			Step<-X.bg$x[2]- X.bg$x[1]
			sum2=sum2+sum(pnbinom(i-para[K], para[3],,para[2])*X.bg$y)*Step 
		}
	}
	if (method=="Poisson")
	{
		sum2=sum2+ sum(ppois(c((cut+1-para[K]):(T-para[K])), para[3])*dpois(c((cut+1):T), para[2]))
	}
	if (method=="NB")
	{
		sum2=sum2+ sum(pnbinom(c((cut+1-para[K]):(T-para[K])), para[3],, para[2])*dnbinom(c((cut+1):T), para[5],,para[4]))
	}
	IPE=round(1-sum2, digits=4)
	return(IPE)
}
