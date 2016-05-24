"ivphiBB6" <-
function(theta,delta,t)
{
	ivphiBB6<-(-log(-((-(exp(-(-log(t))^(1/delta))-1))^(theta)-1)))^(delta)
	resultivphiBB6<-matrix(c(theta,delta,t,ivphiBB6),nrow=1)
}

