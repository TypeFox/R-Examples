rdplot <- function(x, res, f=0.8){
	plot(x, res)
	lines(lowess(x, res,f=f),col=rgb(0.5,0.5,0.5,0.5),lwd=2)
}

