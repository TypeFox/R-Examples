tango.stat<-function(data, listw, zero.policy=FALSE)
{
	p<-scale(data$Expected, center=FALSE, scale=sum(data$Expected))
	r<-scale(data$Observed, center=FALSE, scale=sum(data$Observed))

	rp<-r-p
	T<-lag.listw(listw, rp, zero.policy = zero.policy)
	T<-t(rp)%*%T
	return(T)
}
