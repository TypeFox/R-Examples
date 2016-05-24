F.start.limits <- function( like, expan, w.lo, w.hi, dist ){
#
#   Establish starting value for parameters, and limits passed to the optimizer
#

#   Number of parameters
np <- expan + 1 + 1*(like %in% c("hazrate","uniform"))

w <- w.hi - w.lo


#   No starting values given
if( like == "hazrate" ){
    start <- c(.5*w, 1,rep(0, np - 2))
    low   <- c(0, .01, rep(-Inf, np - 2 ))
    high  <- c(Inf, Inf, rep( Inf, np - 2 ))
    nms <- c("Sigma", "Beta")
    if(expan > 0) nms <- c(nms, paste( "a", 1:(np-2), sep=""))

} else if( like == "halfnorm" ){
    start <- c(sqrt(sum( (dist - w.lo)^2 )/length(dist)), rep(0, np - 1))
    low   <- c(0, rep(-Inf, np - 1 ))
    high  <- c(Inf, rep( Inf, np - 1 ))
    nms <- c("Sigma")
    if(expan > 0) nms <- c(nms, paste( "a", 1:(np-1), sep=""))

} else if( like == "uniform" ){
    start <- c(.1*w, 1, rep(0, np - 2))
    low   <- c(1e-6, 0, rep(-Inf, np - 2 ))
    high  <- c(w, Inf, rep( Inf, np - 2 ))
    nms <- c("Threshold", "Knee")
    if(expan > 0) nms <- c(nms, paste( "a", 1:(np-2), sep=""))

} else if( like == "negexp" ){
    start <- c(1, rep(0, np - 1))
    low   <- c(0, rep(-Inf, np - 1 ))
    high  <- c(Inf, rep( Inf, np - 1 ))
    nms <- c("Beta")
    if(expan > 0) nms <- c(nms, paste( "a", 1:(np-1), sep=""))


} else if( like == "Gamma" ){
    d <- dist[ w.lo <= dist & dist <= w.hi ]
    d <- d[ d > 0 ] # even though 0 is legit, can't take log of it
    s <- log( mean(d, na.rm=TRUE) ) - mean( log(d), na.rm=TRUE )
    s2 <- (s-3)^2 + 24*s
    if( s2 < 0 ) s2 <- 0
    r <- (3 - s + sqrt( s2 )) / (12*s)
    if( r <= 1 ) r <- 1.01
    b <- ( (r-1) / exp(1) )^(r-1) / gamma(r)
    lam <- mean(d,na.rm=TRUE) / (r * b)

    start <- c(r, lam)
    low   <- c(1.0001, 0)
    high  <- c(Inf, Inf)
    nms <- c("Shape", "Scale")


} else {
    #   Assume this is a user-defined likelihood
    fn <- match.fun( paste(like, ".start.limits", sep="") )
    ans <- fn(dist, expan, w.lo, w.hi)
    start <- ans$start
    low <- ans$lowlimit
    high <- ans$highlimit
    nms <- ans$names
}

names(start) <- nms
names(low) <- nms
names(high) <- nms

list( start=start, lowlimit=low, uplimit=high, names=nms )

}
