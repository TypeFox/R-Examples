pgumbel <-
function(q,scale=1,location=0,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- exp(-exp(-(q-location)/scale))
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}

