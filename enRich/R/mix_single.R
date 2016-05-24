mix_single <-function(data1, method=NULL, initialpara=NULL, fixoffset=FALSE, fixk=3,  krange=c(0:10), stopdiff=0.0001)
{
## data1: vector of counts of the data
## method= "Poisson" or "NB"
	
	if (is.null(initialpara)==1)
	{
		if (method=="Poisson")
			initialpara=c(0.1, 10, 1)
		if (method=="NB")
			initialpara=c(0.1, 10, 1,1,1)
	}
	if (fixoffset==FALSE)
	{
      	loglike=rep(-Inf, length(krange))
		para2<-matrix(0, length(krange), (length(initialpara)+2))
      	for (j in 1:length(krange))
      	{
			para1<-mix_offset(data1, method=method, initialpara=initialpara,  k=krange[j], stopdiff=stopdiff)
			if (method=="Poisson")
			{
				para=para1[1:3]
				loglike[j]=mloglike_offsetpois(data1, para1)
            	}
			if (method=="NB")
        		{
				para=para1[1:5]
				loglike[j]=mloglike_offsetNB(data1, para1)
			}
			para2[j,]=c(para1, loglike[j])
		}
		maxj=which.max(loglike)
		maxpara=para2[maxj, ]
		if (method=="Poisson")
		{
			names(maxpara)=c("p", "lambda_S", "lambda_B", "offset k", "loglik")
			colnames(para2)=c("p", "lambda_S", "lambda_B", "offset k", "loglik")
		}
		if (method=="NB")
		{
			names(maxpara)=c("p", "mu_S", "phi_S", "mu_B", "phi_B", "offset k", "loglik")
			colnames(para2)=c("p", "mu_S", "phi_S", "mu_B", "phi_B", "offset k", "loglik")
		}
		results=maxpara
	}
	if (fixoffset==TRUE)
	{
     		para1<-mix_offset(data1, method=method, initialpara=initialpara,  k=fixk, stopdiff=stopdiff)
		if (method=="Poisson")
		{
			loglike1=mloglike_offsetpois(data1, para1)
     		}
		if (method=="NB")
      		{
			loglike1=mloglike_offsetNB(data1, para1)
		}
		maxpara=c(para1, loglike1)
		if (method=="Poisson")
		{
			names(maxpara)=c("p", "lambda_S", "lambda_B", "offset k", "loglik")
		}
		if (method=="NB")
		{
			names(maxpara)=c("p", "mu_S", "phi_S", "mu_B", "phi_B", "offset k", "loglik")
		}
		results=maxpara
	}
	return(results)
}
