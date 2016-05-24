whittermore.stat<-function(data, listw, zero.policy=FALSE)
{
	n<-length(data$Observed)
	r<-scale(data$Observed, center=FALSE, scale=sum(data$Observed) )

	T<-lag.listw(listw, r, zero.policy = zero.policy)
	T<-(n/(n-1))*(t(r)%*%T)

	return(T)
}
