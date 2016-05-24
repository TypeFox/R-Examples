stone.stat<-function(data, region, sorted=FALSE, lambda=NULL)
{
	if(is.null(lambda))
		lambda<-sum(data$Observed)/sum(data$Expected)

	if(!sorted)
	{
		xd<-data$x-data$x[region]
		yd<-data$y-data$y[region]
		dist<-xd*xd+yd*yd

		idx<-order(dist)
	}
	else
	{
		idx<-1:length(data$Observed)
	}
	
	ccoc<-cumsum(data$Observed[idx])/(lambda*cumsum(data$Expected[idx]))
	
	t<-max(ccoc)

	return( c(t, region=which(ccoc==t)) )
}
