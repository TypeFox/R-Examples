"psiGumbel" <-
function(delta,s)
{
		psiGumbel<-exp(-s^(1/delta))
	resultpsiGumbel<-matrix(c(delta,s,psiGumbel),nrow=1)
}

