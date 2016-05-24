dlnorm3 <-
function(x,shape=1,scale=1,thres=0,log=FALSE)
{
	fx <- dlnorm(x-thres,scale,shape)
	if(log) return(log(fx))
	else return(fx)
}

