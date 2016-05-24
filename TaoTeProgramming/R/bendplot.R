"bendplot" <-
function (num=2e4, xdelta=.1, ydelta=.2, sd=1, lwd=2, 
          color="black", seed=NULL) 
{
	if(length(seed)) set.seed(seed)

	plot(bend(num=num, xdelta=xdelta, ydelta=ydelta, sd=sd), type="l",
		axes=FALSE, xlab="", ylab="", lwd=lwd, col=color[1])
}

