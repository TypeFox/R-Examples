"ivphiBB2" <-
function(theta,delta,t)
{
	ivphiBB2<-t^(-theta)-1
	resultivphiBB2<-matrix(c(theta,delta,t,ivphiBB2),nrow=1)
}

