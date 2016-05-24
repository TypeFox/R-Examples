dllog <-
function(x,shape=1,scale=1,log=FALSE)
{
	fx <- dlogis(log(x),location=scale,scale=shape,log=FALSE)/x
	if(log) return(log(fx))
	else return(fx)
}

