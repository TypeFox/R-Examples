hist.plot <-
function(x,...)
{
	hist(x$x.info[,'x'],freq=FALSE,col='gray',main='',xlab=x$data.name,border='white',cex.axis=.8)
	fun <- function(a) distr(a,x$dist,x$par.hat)
	rug(x$x.info[,'x'],col='red',lwd=2)
	curve(fun,add=TRUE,col='red')
}

