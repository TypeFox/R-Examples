pottwhitt.stat <-function(data)
{
	T<-sum(data$Expected)*sum(data$Observed*(data$Observed-1)/data$Expected)

	allO<-sum(data$Observed)
	
	asintmean<-allO*(allO-1)
	asintvar<-2*(length(data$Observed)-1)*asintmean
	pvalue<-pnorm(T, asintmean, sqrt(asintvar))
	pvalue<-min(pvalue, 1-pvalue)
	
	return(list(T=T, asintmean=asintmean, asintvat=asintvar, pvalue=pvalue))
}
