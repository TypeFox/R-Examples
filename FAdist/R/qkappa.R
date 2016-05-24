qkappa <-
function(p,shape=1,scale=1,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- (((scale*p)^(-shape)-1/scale^shape)/shape)^(-1/shape)
	return(xF)
}

