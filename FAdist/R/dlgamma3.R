dlgamma3 <-
function(x,shape=1,scale=1,thres=1,log=FALSE)
{
	fx <- dgamma3(log(x),shape,1/scale,thres)/x
	if(log) return(log(fx))
	else return(fx)
}

