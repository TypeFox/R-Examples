"phiBB6" <-
function(theta,delta,s)
{
	phiBB6<-exp(-(-log(1-(1-exp(-s^(1/delta)))^(1/theta)))^(delta))
	resultphiBB6<-matrix(c(theta,delta,s,phiBB6),nrow=1)
}

