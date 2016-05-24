rlnorm3 <-
function(n,shape=1,scale=1,thres=0)
	rlnorm(n,scale,shape)+thres

