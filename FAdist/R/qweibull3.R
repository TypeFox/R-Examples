qweibull3 <-
function(p,shape,scale=1,thres=0,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- thres+qweibull(p,shape,scale)
	return(xF)
}

