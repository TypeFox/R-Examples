pweibull3 <-
function(q,shape,scale=1,thres=0,lower.tail=TRUE,log.p=FALSE)
	pweibull(q-thres,shape,scale,lower.tail,log.p)

