plot.cband <- function(x,cbands=TRUE,pbands=TRUE,
                       xlab=NULL,ylab=NULL,main=NULL,...) {
#
# Function plot.cband.  To plot data, fitted lines, confidence
# and prediction bands for a mixture of regressions.
# The argument ``x'' is as returned by the function cband().
#

y     <- x$y
xf    <- x$xf
bnds  <- x$bnds
type  <- x$type
theta <- x$theta
int   <- x$intercept
alpha <- x$alpha
x     <- x$x
K     <- length(theta)
ylim  <- range(c(y,unlist(bnds)),na.rm=TRUE)
xt    <- range(x)

i1 <- if(type=='both') 1:2 else 1
i2 <- if(type=='both') 3:4 else 2
if(is.null(main)) {
	tit <- if(type=='both')
		paste('Prediction and confidence bands, level = ',
                               100*(1-alpha),'%.',sep='')
	else if(type=='upper')
		paste('Upper prediction and confidence bands, level = ',
                              100*(1-alpha),'%.',sep='')
	else
		paste('Lower prediction and confidence bands, level = ',
            100*(1-alpha),'%.',sep='')
} else tit <- main
if(is.null(xlab)) xlab <- 'x'
if(is.null(ylab)) ylab <- 'y'

plot(0,0,type='n',xlim=range(x),ylim=ylim,xlab=xlab,ylab=ylab)
points(x,y)
for(j in 1:K) {
	beta <- theta[[j]]$beta
	yt  <- if(int) cbind(1,xt)%*%beta else beta*xt
	lines(xt,yt)
	if(cbands)
		for(i in i1) lines(xf,bnds[[j]][,i],lty=2)
	if(pbands)
		for(i in i2) lines(xf,bnds[[j]][,i],lty=3)
}
title(main=tit)

invisible()
}
