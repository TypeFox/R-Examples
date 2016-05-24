"phiBB1" <-
function(theta,delta,s)
{
		phiBB1<-exp(-(1/theta*log(1+s^(1/delta)))^delta)
		resultphiBB1<-matrix(c(theta,delta,s,phiBB1),nrow=1)
}

