qlgamma3 <-
function(p,shape=1,scale=1,thres=1,lower.tail=TRUE,log.p=FALSE)
{
	if(log.p) p <- exp(p)
	if(!lower.tail) p <- 1 - p
	xF <- exp(qgamma3(p,shape,1/scale,thres))
	return(xF)
}

