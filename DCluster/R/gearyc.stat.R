gearyc.stat<-function(data, applyto="SMR", ...)
{
	n<-length(data$Observed)

	if(applyto != "SMR")
	{
		Z<- data$Observed - data$Expected
	}
	else
	{
		Z<- data$Observed/data$Expected 
		Z[!is.finite(Z)]<-0
	}

	return(spdep::geary(x=Z, ...)$C)
}
