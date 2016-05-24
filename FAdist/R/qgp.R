qgp <-
function(p,shape=1,scale=1,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- scale/shape*(1-(1-p)^shape)
	return(xF)
}

