qllog <-
function(p,shape=1,scale=1,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- exp(qlogis(p,location=scale,scale=shape))
	return(xF)
}

