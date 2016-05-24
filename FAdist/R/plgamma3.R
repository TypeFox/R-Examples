plgamma3 <-
function(q,shape=1,scale=1,thres=1,lower.tail=TRUE,log.p=FALSE)
{
	Fx <- pgamma3(log(q),shape,1/scale,thres)
	if(!lower.tail) Fx <- 1 - Fx
	if(log.p) Fx <- log(Fx)
	return(Fx)
}

