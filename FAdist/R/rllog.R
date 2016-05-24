rllog <-
function(n,shape=1,scale=1)
	exp(rlogis(n,location=scale,scale=shape))

