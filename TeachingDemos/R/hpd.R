# These functions were written by Greg Snow (greg.snow@ihc.com)
# They are free to use, but come with no warrenty whatsoever
# use at your own risk (not that I can think of anything bad that
# they would do).


hpd <- function(posterior.icdf, conf=0.95, tol=0.00000001,...){
	conf <- min( conf, 1-conf )
	f <- function(x,posterior.icdf,conf,...){
		posterior.icdf(1-conf+x,...) - posterior.icdf(x,...)
	}
	out <- optimize(f, c(0,conf), posterior.icdf = posterior.icdf,
                        conf=conf, tol=tol, ...)
	return( c( posterior.icdf(out$minimum,...), 
	           posterior.icdf(1-conf+out$minimum,...) ) )
}

emp.hpd <- function(x, conf=0.95){
	conf <- min(conf, 1-conf)
	n <- length(x)
	nn <- round( n*conf )
	x <- sort(x)
	xx <- x[ (n-nn+1):n ] - x[1:nn]
	m <- min(xx)
	nnn <- which(xx==m)[1]
	return( c( x[ nnn ], x[ n-nn+nnn ] ) )
}
