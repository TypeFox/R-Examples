ESW <- function( obj ){
#
#   obj = a dfunc object.  It may optionally contain a g0 component.
#       if no g0, assume g0 = 1
#
#   Note, For the classic distance functions, the intergral from 0 to w of g(x) = ESW 
#   is the ratio g(x) / f(x) for any x.  If we know g(0), we can compute f(0), and 
#   take the ratio g(0) / f(0).  However, in general and for other distance functions
#   we will do numerical integration.  This is a bit slower, and a bit less accurate in 
#   some cases, but is more general.  This allows user defined distance functions to be 
#   added easily. 

like <- match.fun(paste( obj$like.form, ".like", sep=""))

if( (obj$like.form == "hazrate") & (obj$x.scl == obj$w.lo) ){
    x <- seq( obj$w.lo + 1e-6*(obj$w.hi - obj$w.lo), obj$w.hi, length=200)
} else {
    x <- seq( obj$w.lo, obj$w.hi, length=200)
}

y <- like( obj$parameters, x - obj$w.lo, series=obj$series, expansions=obj$expansions, w.lo = obj$w.lo, w.hi=obj$w.hi )
#print(y)


if( is.null( obj$g.x.scl ) ){
    #   Assume g0 = 1
    g.at.x0 <- 1
    x0 <- 0
    warning("g0 unspecified.  Assumed 1.")
} else {
    g.at.x0 <- obj$g.x.scl
    x0 <- obj$x.scl
}
f.at.x0 <- like( obj$parameters, x0 - obj$w.lo, series=obj$series, expansions=obj$expansions, w.lo=obj$w.lo, w.hi=obj$w.hi )

y <- y * g.at.x0 / f.at.x0


esw <- (x[3] - x[2]) * sum(y[-length(y)]+y[-1]) / 2   # Trapazoid rule.  Use x[3] and x[2] because for hazard rate, x[1] is not evenly spaced with rest

### This is correct. It works provided g(x) and x are stored in obj. 
### This routine evaluates f(x) using like(), 
### and then takes ratio of g(x) / f(x). 
#esw2 <- g.at.x0 / f.at.x0

#print( c(esw2, esw))

esw

}
