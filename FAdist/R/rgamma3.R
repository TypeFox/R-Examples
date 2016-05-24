rgamma3 <-
function(n,shape=1,scale=1,thres=0)
	rgamma(n,shape,1/scale)+thres

