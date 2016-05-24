pkappa <-
function(q,shape=1,scale=1,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- q/scale*(shape+(q/scale)^shape)^(-1/shape)
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}

