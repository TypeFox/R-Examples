achisq.stat<-function(data, lambda=NULL)
{

	df<-length(data$Observed)


	#If internal standardization was  used then lambda=1 and df=n-1
	if(sum(data$Observed)==sum(data$Expected))
	{
		lambda<-1
		df<-df-1
	}
	else 
	{
		#If lambda is unknown then we must slightly modify E_i
		if(is.null(lambda))
		{
			lambda<-sum(data$Observed)/sum(data$Expected)
			df<-df-1
		}
	}
		
	Elambda<-data$Expected*lambda

	T<-sum((data$Observed-Elambda)^2/Elambda)
	pvalue<-pchisq(T, df, lower.tail=FALSE)
	
	return( list(T=T, df=df, pvalue=pvalue) )
}
