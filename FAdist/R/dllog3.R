dllog3 <-
function(x,shape=1,scale=1,thres=0,log=FALSE)
{
	fx <- dllog(x-thres,shape,scale)
	if(log) return(log(fx))
	else return(fx)
}

