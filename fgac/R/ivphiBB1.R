"ivphiBB1" <-
function(theta,delta,t)
{
	ivphiBB1<-(exp(theta*(-log(t))^(1/delta))-1)^delta
	resultivphiBB1<-matrix(c(theta,delta,t,ivphiBB1),nrow=1)
}

