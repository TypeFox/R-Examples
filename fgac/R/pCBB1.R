"pCBB1" <-
function(theta,delta,s,t)
{
	pc<-pcopula1(theta,delta,psiGumbel,phiBB1,ivpsiGumbel,ivphiBB1,s,t);
	resu<-pc
}

