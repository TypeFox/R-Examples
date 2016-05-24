pgev <-
function(q,shape=1,scale=1,location=0,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- exp(-(1+shape*((q-location)/scale))^(-1/shape))
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}

