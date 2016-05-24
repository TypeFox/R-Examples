"psiKS" <-
function(delta,s)
{
	psiKS<-(1+s)^(-1/delta)
	resultpsiKS<-matrix(c(delta,s,psiKS),nrow=1)	
}

