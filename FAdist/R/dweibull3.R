dweibull3 <-
function(x,shape,scale=1,thres=0,log=FALSE)
	dweibull(x-thres,shape,scale,log)

