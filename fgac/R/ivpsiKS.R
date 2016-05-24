"ivpsiKS" <-
function(delta,t)
{
	ivpsiKS<-t^(-delta)-1
	resultivpsiKS<-matrix(c(delta,t,ivpsiKS),nrow=1)
}

