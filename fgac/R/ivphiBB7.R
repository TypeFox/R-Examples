"ivphiBB7" <-
function(theta,delta,t)
{
	ivphiBB7<-(-((-((-(log(t)-1))^(-1/delta)-1))^(theta)-1))^(-delta)-1
	resultivphiBB7<-matrix(c(theta,delta,t,ivphiBB7),nrow=1)
}

