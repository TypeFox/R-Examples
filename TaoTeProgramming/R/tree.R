"tree" <-
function (branches=30, seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	canvas()
	# trunk
	bundle(10, orig=c(.5, .08), dest=c(.5, .95), width=c(.08, .1, .02, .1))
	bundle(20, orig=c(.5, .1), dest=c(.5, .75), width=c(.1, .1, .04, .1))
	bundle(20, orig=c(.5, .1), dest=c(.5, .65), width=c(.1, .1, .05, .2))

	for(i in 1:branches) {
		height <- runif(1, .5, .95)
		extent <- (1 - height) * sign(runif(1, -1, 1))
		bundle(10, orig=c(.5, height), 
			dest=c(.5 + extent, height - runif(1, 0, .2) * 
			abs(extent)), width=c(.01, .1, .05, .05))
	}

}
