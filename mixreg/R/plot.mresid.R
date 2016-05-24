plot.mresid <- function(x,vs.fit=FALSE,whichx=1,shape='disc',
                        ngon=20,size=1,xlab=NULL,...) {
	K      <- ncol(x$resid)
	resid  <- c(x$resid)
	gamma  <- c(x$gamma)
	if(vs.fit) {
		x <- x$y - resid
		if(is.null(xlab)) xlab <- 'fitted values'
	}
	else {
		x <- rep(as.matrix(x$x)[,whichx],K)
		if(is.null(xlab)) xlab <- 'x'
	}
	ishape <- pmatch(shape,c('disc','diamond','square'))
	if(is.na(ishape))
		stop(paste('Shape ',shape,' not recognized.',sep=''))
	plot(x,resid,type='n',xlab=xlab,ylab='residuals')
        uin <- par.uin()

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
                                     sym$y*gamma[i],...)
	invisible()
}
