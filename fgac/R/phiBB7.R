"phiBB7" <-
function(theta,delta,s)
{
	phiBB7<-exp(1-(1-(1-(1+s)^(-1/delta))^(1/theta))^(-delta))
	resultphiBB7<-matrix(c(theta,delta,s,phiBB7),nrow=1)
}

