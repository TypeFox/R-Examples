qq.mix <- function(object,xlim=NULL,ylim=NULL,shape='disc',ngon=20,size=1,...)
{
	K      <- ncol(object$resid)
	resid  <- c(object$resid)
	gamma  <- c(object$gamma)
	ind    <- order(resid)
	resid  <- resid[ind]
	gamma  <- gamma[ind]
	N      <- length(resid)
	x      <- qnorm((0.5+(1:N))/(N+1))

	ishape <- pmatch(shape,c('disc','diamond','square'))
	if(is.na(ishape))
		stop(paste('Shape ',shape,' not recognized.',sep=''))
	if(is.null(xlim)) xlim <- range(x)
	if(is.null(ylim)) ylim <- range(resid)
	plot(x,resid,type='n',xlim=xlim,ylim=ylim,
             xlab='normal quantiles',ylab='empirical quantiles')
	uin  <- par.uin() # user units per inch
	if(ishape==1) {
		phi <- c(seq(0,2*pi,length=ngon),0)
		scale <- 0.03*size/uin
		sym <- list(x=cos(phi)*scale[1],y=sin(phi)*scale[2])
	} else if(ishape==2) {
		scale <- 0.03*size/uin
		sym <- list(x=c(1,0,-1,0,1)*scale[1],y=c(0,1,0,-1,0)*scale[2])
	} else {
		scale <- 0.03*size/(sqrt(2)*uin)
		sym <- list(x=c(1,-1,-1,1,1)*scale[1],y=c(1,1,-1,-1,1)*scale[2])
	}
	N <- length(x)
	for(i in 1:N) polygon(x[i] + sym$x*gamma[i],resid[i] +
                                     sym$y*gamma[i],err=-1,...)
	invisible()
}
