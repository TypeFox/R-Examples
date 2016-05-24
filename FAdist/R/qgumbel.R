qgumbel <-
function(p,scale=1,location=0,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- location-scale*log(-log(p))
	return(xF)
}

