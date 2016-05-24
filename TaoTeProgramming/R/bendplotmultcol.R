"bendplotmultcol" <-
function (num=2e4, xdelta=.1, ydelta=.2, sd=1, lwd=2, color=colors(),
	cnum=100, seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	bobj <- bend(num=num, xdelta=xdelta, ydelta=ydelta, sd=sd)
	plot(bobj, type="n", axes=FALSE, xlab="", ylab="", lwd=lwd)
	block <- seq(1, floor(num / cnum))
	offset <- max(block) * 0:cnum
	for(i in 1:cnum) {
		lines(bobj[offset[i] + block, ], col=safesample(color))
	}
}

