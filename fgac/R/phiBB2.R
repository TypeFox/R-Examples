"phiBB2" <-
function(theta,delta,s)
{
	phiBB2<-(1+s)^(-1/theta)
	resultphiBB2<-matrix(c(theta,delta,s,phiBB2),nrow=1)
}

