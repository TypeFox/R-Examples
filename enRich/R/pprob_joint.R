pprob_joint <-
function(data1, para, method=NULL)
{
## data1: the counts of two or more experiments
## para: fitting results of these experiments
	N=nrow(data1)
	Nsr=ncol(data1)

	pr1=rep(0, N)
	pr0=rep(0, N)
   	commonp1=0

   	for(j in 1:Nsr)
   	{
		if (method=="Poisson")
		{
			pr1=pr1+ifelse(data1[,j]<para[(j-1)*4+4], log(0.00001), dpois(data1[,j]-para[(j-1)*4+4], para[(j-1)*4+2], log=TRUE))
		   	pr0=pr0+dpois(data1[,j], para[(j-1)*4+3], log=TRUE)
			commonp1=commonp1+para[(j-1)*4+1]
		}
		if (method=="NB")
		{
			pr1=pr1+ifelse(data1[,j]<para[(j-1)*6+6], log(0.00001), dnbinom(data1[,j]-para[(j-1)*6+6], para[(j-1)*6+3],,para[(j-1)*6+2], log=TRUE))
		   	pr0=pr0+dnbinom(data1[,j], para[(j-1)*6+5],,para[(j-1)*6+4], log=TRUE)
			commonp1=commonp1+para[(j-1)*6+1]
		}
   	}
	commonp1=commonp1/Nsr
	pr1=pr1+log(commonp1)
	pr0=pr0+log(1-commonp1)
	ppr1=1/(1+exp(pr0-pr1))
	ppr0=1-ppr1
	pp=list(px1=ppr1, px0=ppr0)
	return(pp)
}
