"phiBB3" <-
function(theta,delta,s)
{
	phiBB3<-exp(1-exp(delta*(theta^(-1)*log(1+s))^(1/delta)))
	resultphiBB3<-matrix(c(theta,delta,s,phiBB3),nrow=1)
}

