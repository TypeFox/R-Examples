plot.dfunc <- function( x, include.zero=FALSE, nbins="Sturges", ... ){
#
#   Plot method for distance functions.
#
#   input: x = object of class "dfunc"
#   include.zero = whether or not to plot distance function at 0.

#   changed the number of plotting points to 200 - jg
x.dist <- x$dist[ (x$w.lo <= x$dist) & (x$dist <= x$w.hi) ]
cnts <- hist( x.dist, plot=FALSE, breaks=nbins )
xscl <- cnts$mid[2] - cnts$mid[1]

#   Gotta add bars on the left if first bar is not at w.lo.  I.e., if first 
#   bar is zero.  Zero bars at top end are not a problem, but low end are because
#   barplot just plots bars, not coordinates
if( cnts$breaks[1] > x$w.lo ){
    # do the hist again, this time specifying breaks exactly
    brks <- seq(x$w.lo, x$w.hi, by=xscl)
    brks <- c(brks, brks[length(brks)] + xscl )   # make sure last bin goes outside range of data
    cnts <- hist( x.dist, plot=FALSE, breaks=brks, include.lowest=TRUE )
}


like <- match.fun( paste( x$like.form, ".like", sep=""))

x.seq <- seq( x$w.lo, x$w.hi, length=200)

y <- like( x$parameters, x.seq - x$w.lo, series=x$series, expansions=x$expansions, w.lo=x$w.lo, w.hi=x$w.hi )

if( include.zero & x$like.form == "hazrate" ){
    x.seq[1] <- x$w.lo
}

if( is.null( x$g.x.scl ) ){
    #   Assume g0 = 1
    g.at.x0 <- 1
    x0 <- 0
    warning("g0 unspecified.  Assumed 1.")
} else {
    g.at.x0 <- x$g.x.scl
    x0 <- x$x.scl
}
f.at.x0 <- like( x$parameters, x0 - x$w.lo, series=x$series, expansions=x$expansions, w.lo=x$w.lo, w.hi=x$w.hi )
if(is.na(f.at.x0) | (f.at.x0 <= 0)){
    #   can happen when parameters at the border of parameter space
    yscl <- 1.0
    warning("Y intercept missing or zero. One or more parameters likely at their boundaries. Caution.")
} else {
    yscl <- g.at.x0 / f.at.x0
}


ybarhgts <- cnts$density * yscl
y <- y * yscl

y.finite <- y[ y < Inf ]
y.lims <- c(0, max( g.at.x0, ybarhgts, y.finite, na.rm=TRUE ))

if( include.zero ){
    x.limits <- c(0 , max(x.seq))
} else {
    x.limits <- range(x.seq) 
}


bar.mids <- barplot( ybarhgts, width=xscl, space=0, density=0, ylim=y.lims, xlim=x.limits, border="blue", ... )   # real x coords are 1, 2, ..., nbars
xticks <- axTicks(1)
axis( 1, at=xticks,  labels=xticks, line=.5 )
title( xlab="Distance", ylab="Probability of detection" )
if( !("main" %in% names(c(...))) ){
    # Put up a default title containing liklihood description
    if( x$expansions == 0 ){
        title(main=paste( x$like.form, ", ", x$expansions, " expansions", sep=""))
    } else {
        title(main=paste( x$like.form, ", ", x$series, " expansion, ", x$expansions, " expansions", sep=""))
    }
}

#   These 3 lines plot a polygon for the density function
#x.poly <- c(0, x, x[length(x)] )
#y.poly <- c(0, y, 0)
#polygon( x.poly, y.poly, density=15, border="red", lwd=2 )

#   This places a single line over the histogram
lines( x.seq, y, col="red", lwd=2 )

#   These two add vertical lines at 0 and w
lines( rep(x.seq[1], 2), c(0,y[1]), col="red", lwd=2 )
lines( rep(x.seq[length(x.seq)], 2), c(0,y[length(x.seq)]), col="red", lwd=2 )

#   print area under the curve
area <- ESW( x )
#area2 <- (x[3] - x[2]) * sum(y[-length(y)]+y[-1]) / 2   # use x[3] and x[2] because for hazard rate, x[1] is not evenly spaced with rest
#print(c(area,area2))
text( max(x.seq), max(y.lims)-0.025*diff(y.lims), paste("ESW =", round(area,3)), adj=1)

#   If model did not converge, print a message on the graph.
if( x$convergence != 0 ){
    if( x$convergence == -1 ){
        mess <- "Solution failure"
    } else {
        mess <- "Convergence failure"
    }
    text( mean(x.seq), mean(y.lims), mess, cex=3, adj=.5, col="red")
    text( mean(x.seq), mean(y.lims), paste("\n\n\n", x$fit$message, sep=""), cex=1, adj=.5, col="black")
} else if( is.na(area) | (area > x$w.hi) ){
    #   invalid scaling, g0 is wrong
    mess <- "Scaling failure"
    text( mean(x.seq), mean(y.lims), mess, cex=3, adj=.5, col="red")
    mess <- paste("Check g0=", round(g.at.x0,2), "assumption")
    text( mean(x.seq), mean(y.lims), paste("\n\n\n", mess,sep=""), cex=1, adj=.5, col="black")
}
#x$xscl.plot <- xscl   # gonna need this to plot something on the graph.
x$yscl <- yscl   # this is g(x) / f(x).  Might want this later.

invisible(x)
}
