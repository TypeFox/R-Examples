plot.zahb <- function(x, add=FALSE, ...) {
    plotAstro(x$data$logTe, x$data$logL, yi=0.1, xi=0.1, xlab=expression(log ~ italic(T)[eff]), ylab=expression(log ~ italic(L/L)[sun]), revX=TRUE, xmt=4, ymt=4, add=add, ...)
}

plot.trk <- function(x, add=FALSE, ...) {
    plotAstro(x$data$logTe, x$data$logL, yi=0.5, xi=0.05, xlab=expression(log ~ italic(T)[eff]), ylab=expression(log ~ italic(L/L)[sun]), revX=TRUE, xmt=4, ymt=4, add=add, ...)
}

plot.hb <- function(x, add=FALSE, ...) {
   sel <- x$data$time >= 6
   plotAstro(x$data$logTe[sel], x$data$logL[sel], yi=0.2, xi=0.05, xlab=expression(log ~ italic(T)[eff]), ylab=expression(log ~ italic(L/L)[sun]), revX=TRUE, xmt=4, ymt=4, add=add, ...)
}

plot.iso <- function(x, add=FALSE, ...) {
    plotAstro(x$data$logTe, x$data$logL, yi=1, xi=0.05, xlab=expression(log ~ italic(T)[eff]), ylab=expression(log ~ italic(L/L)[sun]), revX=TRUE, xmt=4, ymt=4, add=add, ...)
}

################################################

plot.trkset <- function(x, add=FALSE, col=1, lty=1, xlim=NULL, ylim=NULL, ...) {
    n <- length(x)
    col <- rep(col, length.out=n)
    lty <- rep(lty, length.out=n)

    if(is.null(xlim)) {
        rgTe <- range(sapply(x, function(y) range(y$data$logTe)))

    } else {
        rgTe <- xlim
    }

    if(is.null(ylim)) {
        rgL <- range(sapply(x, function(y) range(y$data$logL)))
    } else {
        rgL <- ylim
    }

    xi <- (rgTe[2]-rgTe[1])/5
    yi <- (rgL[2]-rgL[1])/5
    
    plotAstro(x[[1]]$data$logTe, x[[1]]$data$logL, yi=yi, xi=xi, xlab=expression(log ~ italic(T)[eff]), ylab=expression(log ~ italic(L/L)[sun]), xlim=rgTe, ylim=rgL, revX=TRUE, xmt=4, ymt=4, add=add, col=col[1], lty=lty[1], ...)
    
    for(i in 2:n) 
        plotAstro(x[[i]]$data$logTe, x[[i]]$data$logL, add=TRUE, col=col[i], lty=lty[i], ...)
}

plot.isoset <- function(x, add=FALSE, col=1, lty=1, xlim=NULL, ylim=NULL, ...) {
    n <- length(x)
    col <- rep(col, length.out=n)
    lty <- rep(lty, length.out=n)

    if(is.null(xlim)) {
        rgTe <- range(sapply(x, function(y) range(y$data$logTe)))

    } else {
        rgTe <- xlim
    }

    if(is.null(ylim)) {
        rgL <- range(sapply(x, function(y) range(y$data$logL)))
    } else {
        rgL <- ylim
    }

    xi <- (rgTe[2]-rgTe[1])/5
    yi <- (rgL[2]-rgL[1])/5
    
    plotAstro(x[[1]]$data$logTe, x[[1]]$data$logL, yi=yi, xi=xi, xlab=expression(log ~ italic(T)[eff]), ylab=expression(log ~ italic(L/L)[sun]), xlim=rgTe, ylim=rgL, revX=TRUE, xmt=4, ymt=4, add=add, col=col[1], lty=lty[1], ...)

    for(i in 2:n) 
        plotAstro(x[[i]]$data$logTe, x[[i]]$data$logL, add=TRUE, col=col[i], lty=lty[i], ...)
}


plot.hbset <- function(x, add=FALSE, col=1, lty=1, xlim=NULL, ylim=NULL, ...) {
    n <- length(x)
    col <- rep(col, length.out=n)
    lty <- rep(lty, length.out=n)

    if(is.null(xlim)) {
        rgTe <- range(sapply(x, function(y) range(y$data$logTe)))

    } else {
        rgTe <- xlim
    }

    if(is.null(ylim)) {
        rgL <- range(sapply(x, function(y) range(y$data$logL)))
    } else {
        rgL <- ylim
    }
    
    xi <- (rgTe[2]-rgTe[1])/5
    yi <- (rgL[2]-rgL[1])/5
    
    sel <- x[[1]]$data$time >= 6
    plotAstro(x[[1]]$data$logTe[sel], x[[1]]$data$logL[sel], yi=yi, xi=xi, xlab=expression(log ~ italic(T)[eff]), ylab=expression(log ~ italic(L/L)[sun]), xlim=rgTe, ylim=rgL, revX=TRUE, xmt=4, ymt=4, add=add, col=col[1], lty=lty[1], ...)
    
    for(i in 2:n) {
        sel <- x[[i]]$data$time >= 6
        plotAstro(x[[i]]$data$logTe[sel], x[[i]]$data$logL[sel], add=TRUE, col=col[i], lty=lty[i], ...)
  }
}

###############################################

plotAstro <- function(x, y, type="l", xlab="", ylab="", xi=(max(x)-min(x))/5, yi=(max(y)-min(y))/5, xmt=3, ymt=3, revX=FALSE, revY=FALSE, xlim=NULL, ylim=NULL, cex=1.0, cex.axis=1.3, cex.lab=1.5, add=FALSE, ...)
{
  # minor ticks length
    TCL <- par("tcl")/2
    par(mar=c(4.5,4.7,1,1) + 0.1)
    if(is.null(xlim)) {
        X <- range(x)  }
    else {
        X <- xlim
    }
    if(is.null(ylim)) {
        Y <- range(y)
    }
    else {
        Y <- ylim
    }
  
# number of major ticks
    nx <- ceiling((max(X)-min(X))/xi)
    ny <- ceiling((max(Y)-min(Y))/yi)

# major ticks location 
    xt <- pretty(X, nx)
    yt <- pretty(Y, ny)

# minor ticks location
    stepx <- (xt[2]-xt[1])/(xmt+1)
    stepy <- (yt[2]-yt[1])/(ymt+1)
    xtm <- seq(min(xt), max(xt), by=stepx)
    ytm <- seq(min(yt), max(yt), by=stepy)  
    xtm <- xtm[ ! xtm %in% xt ]
    ytm <- ytm[ ! ytm %in% yt ]

# X and Y ranges  
    xr <- c(min(X,xt), max(X,xt))
    yr <- c(min(Y,yt), max(Y,yt))
    if(revX) {
        xr <- rev(xr)
        xt <- rev(xt)
        xtm <- rev(xtm)
    }
    if(revY) {
        yr <- rev(yr)
        yt <- rev(yt)
        ytm <- rev(ytm)
    }
  
    if(!add) {
# plot without axes
        plot(x,y, type=type, xlab=xlab, ylab=ylab, xlim=xr, ylim=yr, axes=FALSE, cex.lab=cex.lab, cex=cex, ...)
        box()
        
 # add axes
        axis(1, at=xt, lwd=0, lwd.ticks=1, cex.axis=cex.axis)
        axis(1, at=xtm, labels=NA, lwd=0, lwd.ticks=1, tcl=TCL)
        axis(2, at=yt, lwd=0, lwd.ticks=1, cex.axis=cex.axis)
        axis(2, at=ytm, labels=NA, lwd=0, lwd.ticks=1, tcl=TCL)
        
        axis(3, at=xt, labels=NA, lwd=0, lwd.ticks=1)
        axis(3, at=xtm, labels=NA, lwd=0, lwd.ticks=1, tcl=TCL)
        axis(4, at=yt, labels=NA, lwd=0, lwd.ticks=1)
        axis(4, at=ytm, labels=NA, lwd=0, lwd.ticks=1, tcl=TCL)
    } else {
        if(type != "l")
            points(x,y, cex=cex, ...)
        else
            lines(x,y, ...)
    }
}
