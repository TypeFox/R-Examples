pgp <-
function(q,shape=1,scale=1,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- 1-(1-shape*q/scale)^(1/shape)
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}

