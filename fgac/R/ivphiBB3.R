"ivphiBB3" <-
function(theta,delta,t)
{
	ivphiBB3<-exp(theta*(delta^(-1)*log(-(log(t)-1)))^(delta))-1
	resultivphiBB3<-matrix(c(theta,delta,t,ivphiBB3),nrow=1)
}

