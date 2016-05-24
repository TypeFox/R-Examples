plot.SPIn <-
function(x, ...){
	xx <- x$x
	h <- hist(xx,freq=F)
	dens <- x$dens
	hist(xx,freq=F,ylim=c(-.05,max(h$density,dens$y)),main='Histogram')
	lines(dens)
	l <- x$spin[1]
	u <- x$spin[2]
	lines(c(l,u),rep(-.025,2),col='red')
	lines(rep(l,2),c(-.028,-.022),col='red')
#	lines(c(l,l+.05),rep(-.022,2),col='red')
#	lines(c(l,l+.05),rep(-.028,2),col='red')
	lines(rep(u,2),c(-.028,-.022),col='red')
#	lines(c(u,u-.05),rep(-.022,2),col='red')
#	lines(c(u,u-.05),rep(-.028,2),col='red')
	d <- (max(h$breaks)-min(h$breaks))*.04
	text(l-d,-.025,round(l,2))
	text(u+d,-.025,round(u,2))
	conf <- (1 - x$conf)/2
	l <- quantile(xx,conf)
	u <- quantile(xx,1-conf)
	lines(c(l,u),rep(-.05,2))
	lines(rep(l,2),c(-.053,-.047))
#	lines(c(l,l+.05),rep(-.047,2))
#	lines(c(l,l+.05),rep(-.053,2))
	lines(rep(u,2),c(-.053,-.047))
#	lines(c(u,u-.05),rep(-.047,2))
#	lines(c(u,u-.05),rep(-.053,2))
	text(l-d,-.05,round(l,2))
	text(u+d,-.05,round(u,2))
	legend('topright',c('spin','central'),col=c('red','black'),lty=c(1,1))
}
