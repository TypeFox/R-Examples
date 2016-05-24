qllog3 <-
function(p,shape=1,scale=1,thres=0,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- qllog(p,shape,scale)+thres
	return(xF)
}

