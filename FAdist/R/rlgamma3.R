rlgamma3 <-
function(n,shape=1,scale=1,thres=1)
	exp(rgamma3(n,shape,1/scale,thres))

