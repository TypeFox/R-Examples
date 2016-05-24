dgamma3 <-
function(x,shape=1,scale=1,thres=0,log=FALSE)
{
	fx <- dgamma(x-thres,shape,1/scale)
	if(log) return(log(fx))
	else return(fx)
}

