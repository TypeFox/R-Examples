pllog <-
function(q,shape=1,scale=1,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- plogis(log(q),location=scale,scale=shape)
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}

