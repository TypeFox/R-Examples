pllog3 <-
function(q,shape=1,scale=1,thres=0,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- pllog(q-thres,shape,scale)
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}

