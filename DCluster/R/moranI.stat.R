moranI.stat<-function(data, applyto="SMR", ...)
{
	if(applyto != "SMR")
	{
		Z<- data$Observed - data$Expected
	}
	else
	{
		Z<- data$Observed/data$Expected 
		Z[!is.finite(Z)]<-0
	}

	return(spdep::moran(x=Z,...)$I)
}
