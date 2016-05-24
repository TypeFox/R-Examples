pprob <-function(data1, para, method=NULL)
{
## data1: the counts of the experiment
## para: the fitting results of that experiment
## method: "Poisson" or "NB"

	#1. Check data, method and initialpara are consistent
	if (sum(method=="Poisson")+sum(method=="NB")==0)
	{
		stop('method must be given, either "Poisson" or "NB"!', call.=FALSE)
	}
	if (method=="Poisson")
	{
		pr1=log(para[1])+ifelse(data1<para[4], log(0.00001), dpois(data1-para[4], para[2], log=TRUE))
     		pr0=log(1-para[1])+dpois(data1, para[3], log=TRUE)
	}
	if (method=="NB")
	{
		pr1=log(para[1])+ifelse(data1<para[6], log(0.00001), dnbinom(data1-para[6], para[3], ,para[2], log=TRUE))
     		pr0=log(1-para[1])+dnbinom(data1, para[5],,para[4], log=TRUE)
	}
	ppr1=1/(1+exp(pr0-pr1))
     	ppr0=1-ppr1
	pp=list(px1=ppr1, px0=ppr0)
     	return(pp)
}
