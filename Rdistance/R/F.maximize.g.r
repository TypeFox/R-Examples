F.maximize.g <- function( fit ){
#
#   Maximize the distance function in fit.  That is, find x such that g(x) is at its
#   maximum.  G is smooth, so this is easy for nlminb.
#

g.neg <-  function(x, params, like, w.lo=0, w.hi=max(dist), series, expansions=0){

    
    f.like <- match.fun(paste( like, ".like", sep=""))

    g.x <- f.like( params, x, w.lo=w.lo, w.hi=w.hi, series, expansions )

    -g.x * 100000
}

x.start <- (fit$w.lo + fit$w.hi) / 2

x.max <- nlminb(x.start, g.neg,  params = fit$parameters, w.lo=fit$w.lo, w.hi=fit$w.hi, like=fit$like.form,
    expansions=fit$expansions, series=fit$series, lower=fit$w.lo, upper=fit$w.hi)


if( x.max$convergence != 0 ){
    warning(paste("Maximum of g() could not be found. Message=", x.max$message))
    x.max = NA
} else {
    x.max <- x.max$par
}

x.max

}
