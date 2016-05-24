
denstrip <- function(x, dens, at, width, horiz=TRUE, colmax, colmin="white", scale=1, gamma=1, 
                     ticks=NULL, tlen=1.5, twd, tcol, mticks=NULL, mlen=1.5, mwd, mcol, 
                     lattice=FALSE,
                     ...)
{
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(dens)) {
        ## x assumed to be a sample from a distribution, density estimated
        de <- density(x, ...)
        x <- de$x; dens <- de$y
    }
    else {
        if (!is.numeric(dens)) stop("\'dens\' must be numeric")
        if (length(dens) != length(x)) stop("Lengths of \'dens\' and \'x\' must be the same")
        dens <- dens[order(x)]
        x <- sort(x)
    }
    if (lattice) {
        rect.fn <- panel.rect; seg.fn <- panel.segments
        default.width <- diff(current.panel.limits()[[if(horiz) "ylim" else "xlim"]]) / 30
        default.colmax <- trellis.par.get("add.line")$col
        default.twd <- trellis.par.get("add.line")$lwd; default.mwd <- trellis.par.get("add.line")$lwd*2
    }
    else {
        rect.fn <- rect; seg.fn <- segments;
        default.width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 30
        default.colmax <- par("fg")
        default.twd <- par("lwd"); default.mwd <- par("lwd")*2
    }
    if (missing(width)) width <- default.width
    if (missing(colmax)) colmax <- default.colmax
    if (missing(twd)) twd <- default.twd
    if (missing(mwd)) mwd <- default.mwd
    if (missing(tcol)) tcol <- colmax
    if (missing(mcol)) mcol <- colmax
    dens <- dens / max(dens) * scale
    n <- length(x)
    rgbmax <- col2rgb(colmax, alpha=TRUE)
    rgbmin <- if (colmin=="transparent") c(col2rgb(colmax, alpha=FALSE), 0) else col2rgb(colmin, alpha=TRUE)
    if (gamma <= 0) stop("gamma must be greater than 0")
    p <- dens[1:(n-1)] ^ gamma
    if (colmin=="transparent") 
        cols <- rgb(p*rgbmax[1] + (1 - p)*rgbmin[1],
                    p*rgbmax[2] + (1 - p)*rgbmin[2],
                    p*rgbmax[3] + (1 - p)*rgbmin[3],
                    alpha=p*rgbmax[4] + (1 - p)*rgbmin[4],
                    maxColorValue=255)
    else 
        cols <- rgb(p*rgbmax[1] + (1 - p)*rgbmin[1],
                    p*rgbmax[2] + (1 - p)*rgbmin[2],
                    p*rgbmax[3] + (1 - p)*rgbmin[3],
                    alpha=rgbmax[4],
                    maxColorValue=255)
    first.col <- c(TRUE, cols[2:(n-1)] != cols[1:(n-2)])
    next.col <- c(first.col,TRUE); next.col[1] <- FALSE    
    if (horiz) {
        xleft <- x[-n][first.col]; xright=x[next.col];
        ybottom <- at-width/2; ytop <- at+width/2
    }
    else {
        xleft <- at-width/2; xright <- at+width/2;
        ybottom <- x[-n][first.col]; ytop <- x[next.col]
    }
    rect.fn(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop,
         border=NA, col = cols[first.col])
    if (!is.null(ticks)){
        if (horiz) { tx0 <- tx1 <- ticks; ty0 <- at-width*tlen/2; ty1 <- at+width*tlen/2 }
        else { tx0 <- at-width*tlen/2; tx1 <- at+width*tlen/2; ty0 <- ty1 <- ticks }
        seg.fn(tx0, ty0, tx1, ty1, lwd=twd, col=tcol)
    }
    if (!is.null(mticks)){
        if (horiz) { tmx0 <- tmx1 <- mticks; tmy0 <- at-width*mlen/2; tmy1 <- at+width*mlen/2  }
        else { tmx0 <- at-width*mlen/2; tmx1 <- at+width*mlen/2; tmy0 <- tmy1 <- mticks }
        seg.fn(tmx0, tmy0, tmx1, tmy1, lwd=mwd, col=mcol)
    }
    invisible()
}

denstrip.normal <- function(mean, sd, log=FALSE, nx=1000, ...){    
    x <-
        if (log) qlnorm(ppoints(nx, a=0), mean, sd)
        else qnorm(ppoints(nx, a=0), mean, sd)
    dens <- if (log) dlnorm(x, mean, sd) else dnorm(x, mean, sd)
    denstrip(x=x, dens=dens, ...)
}

denstrip.legend <- function(x, # central x position 
                            y, # central y position
                            width, len, colmax, colmin="white", gamma=1, horiz=FALSE,
                            max=1, nticks=5, ticks, value.adj=0, cex, main="Density", lattice=FALSE)
{
    if (lattice) {
        poly.fn <- panel.polygon; seg.fn <- panel.segments; text.fn <- panel.text
        default.width <- diff(current.panel.limits()[[if(horiz) "ylim" else "xlim"]]) / 30
        default.len <- diff(trellis.last.object()[[if(horiz) "x.limits" else "y.limits"]]) / 4
        default.colmax <- trellis.par.get("add.line")$col
        default.colmax <- trellis.par.get("background")$col
        default.cex <- trellis.par.get("axis.text")$cex * 0.75
      }
    else {
        poly.fn <- polygon; seg.fn <- segments; text.fn <- text
        default.width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 30
        default.len <- diff(par("usr")[if(horiz) 1:2 else 3:4]) / 4
        default.colmax <- par("fg")
        default.cex <- par("cex") * 0.75
    }
    if (missing(width)) width <- default.width
    if (missing(len)) len <- default.len
    if (missing(colmax)) colmax <- default.colmax
    if (missing(cex)) cex <- default.cex
    if (horiz) {        
        pt <- x; at <- y
        xdim <- len; ydim <- width
    }
    else {
      pt <- y; at <- x
      xdim <- width; ydim <- len
    }
    npoints <- 1000 # number of distinct colors. no need to make this an argument
    dx <- seq(pt - len/2, pt + len/2, length=npoints)
    ddens <- seq(0, max, length=npoints)
    denstrip(x=dx, dens=ddens, at=at, width=width, colmax=colmax, colmin=colmin, gamma=gamma, horiz=horiz, lattice=lattice)    
    poly.fn(x = c(x - xdim/2, x + xdim/2, x + xdim/2, x - xdim/2), # draw box around strip
            y = c(y - ydim/2, y - ydim/2, y + ydim/2, y + ydim/2))
    if (missing(ticks)) {
        ticks.at <- seq(pt - len/2, pt + len/2, length=nticks)
        tick.nos <- signif(seq(0, max, length=nticks), 3)
    }
    else {
        ticks.at <- pt - len/2 + ticks / max * len
        tick.nos <- ticks
    }
    if (horiz) seg.fn(ticks.at, at - width/2, ticks.at, at - width*0.75)
    else seg.fn(at + width/2, ticks.at, at + width*0.75, ticks.at)
    if (horiz) text.fn(ticks.at, at - width - value.adj, tick.nos, cex = cex, pos=1)
    else text.fn(at + width + value.adj, ticks.at, tick.nos, cex = cex, pos=4)
    text.fn(x, y + ydim/2, main, cex = cex, pos=3)
    invisible()
}

densregion <- function(x, ...) UseMethod("densregion")

densregion.default <- function(x, # times we have estimates for (vector)
                               y, # vector of values on y axis to show distinct densities at.  same number of distinct y values for each t.
                               z, # matrix of densities at t times y
                               pointwise=FALSE,
                               nlevels=100, # number of distinct densities to show
                               colmax=par("fg"),
                               colmin="white",
                               scale=1,
                               gamma=1,
                               contour=FALSE,
                               ...
                               )
{
    if (pointwise)
        z <- z/apply(z, 1, max)
    qz <- unique(quantile(z, seq(0,1,length=nlevels+1), type=1))
    nlevels <- length(qz) - 1
    zq <- cut(z, qz, include.lowest=TRUE, labels=seq(length=nlevels))
    zq <- matrix(as.numeric(zq), nrow=nrow(z), ncol=ncol(z))
    dens <- tapply(z, zq, mean) # unique densities to plot 
    rgbmax <- col2rgb(colmax, alpha=TRUE)
    rgbmin <- if (colmin=="transparent") c(col2rgb(colmax, alpha=FALSE), 0) else col2rgb(colmin, alpha=TRUE)
    if (gamma <= 0) stop("gamma must be greater than 0")
    p <- (dens / max(dens)) ^ gamma
    if (colmin=="transparent") 
        cols <- rgb(p*rgbmax[1] + (1 - p)*rgbmin[1],
                    p*rgbmax[2] + (1 - p)*rgbmin[2],
                    p*rgbmax[3] + (1 - p)*rgbmin[3],
                    alpha=p*rgbmax[4] + (1 - p)*rgbmin[4], maxColorValue=255)
    else 
        cols <- rgb(p*rgbmax[1] + (1 - p)*rgbmin[1],
                    p*rgbmax[2] + (1 - p)*rgbmin[2],
                    p*rgbmax[3] + (1 - p)*rgbmin[3],
                    alpha=rgbmax[4],
                    maxColorValue=255)
    z <- z[order(x),order(y)]
    x <- sort(x)
    y <- sort(y)
    .filled.contour(as.double(x), as.double(y), scale*z/max(z), qz/max(z), col = cols)
    if (contour)
      contour(x, y, scale*z/max(z), add=TRUE, ...)
    invisible()
}

densregion.survfit <- function(x, ny=20, ...) {
    if (!inherits(x, "survfit")) stop(deparse(substitute(x)), " should be a survfit object")
    if (is.null(x$conf.type))
        stop("No confidence intervals available in \"", deparse(substitute(x)), "\"")
    identity <- function(x)x # this is in base R from 2.7 
    tr <- switch(x$conf.type, "log-log"=function(x)log(-log(x)), "log" = log, plain=identity)
    invtr <- switch(x$conf.type, "log-log"=function(x)exp(-exp(x)), "log" = exp, plain=identity)
    drop <- is.na(x$surv) | is.na(x$lower)
    surv <- x$surv[!drop]; lower <- x$lower[!drop]; time <- x$time[!drop]
    lsurv <- tr(surv)
    lse <- abs((tr(surv) - tr(lower)) / qnorm(0.975))
    n <- length(time)
    if (n==0) stop("No observed events")
    y <- matrix(nrow=n, ncol=ny) 
    ## "ny" ordinates, based on normal quantiles, where density must be calculated for each event time
    for (i in 1:n)
        y[i,] <- invtr(qnorm(ppoints(ny, a=0), lsurv[i], lse[i]))
    ## but actually calculate density at all y for every time, defining a grid on the whole plot region
    yy <- sort(unique(y))
    z <- matrix(nrow=n, ncol=length(yy))
    for (i in 1:n)
        z[i,] <- dnorm(tr(yy), lsurv[i], lse[i])
    densregion.default(x=time, y=yy, z=z, ...)
    invisible()
}

densregion.normal <- function(x, mean, sd, ny=20, ...)
{
    n <- length(x) 
    if (n != length(mean))
        stop(deparse(substitute(mean)), " and ", deparse(substitute(x)), " should be the same length")
    if (n != length(sd))
        stop(deparse(substitute(sd)), " and ", deparse(substitute(x)), " should be the same length")
    y <- matrix(nrow=n, ncol=ny)
    for (i in 1:n)
        y[i,] <- qnorm(ppoints(ny, a=0), mean[i], sd[i])
    yy <- sort(unique(y))
    z <- matrix(nrow=n, ncol=length(yy))
    for (i in 1:n)
        z[i,] <- dnorm(yy, mean[i], sd[i])    
    densregion.default(x=x, y=yy, z=z, ...)
    invisible()
}

seqToIntervals <- function(x){
    x <- sort(unique(as.integer(x)))
    breaks <- x[c(1, diff(x)) > 1]
    groups <- cut(x, c(min(x)-1, breaks - 0.5, max(x)+1), labels = FALSE)
    ranges <- tapply(x, groups, range)
    res <- do.call("rbind", ranges)
    colnames(res) <- c("from","to")
    return(res)
}

sectioned.density <- function(x, dens, at, width, offset, ny,
                              method=c("kernel","frequency"), nx, horiz=TRUE, up.left = TRUE,
                              colmax, colmin="white", gamma=1, lattice=FALSE, ...)
{
    if (lattice) {
        rect.fn <- panel.rect
        default.width <- diff(current.panel.limits()[[if(horiz) "ylim" else "xlim"]]) / 20
        default.colmax <- trellis.par.get("add.line")$col
    }
    else {
        rect.fn <- rect
        default.width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 20
        default.colmax <- par("fg")
    }
    if (missing(width)) width <- default.width
    if (missing(colmax)) colmax <- default.colmax
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(offset)) offset <- width/3
    if (!up.left) offset <- - offset
    if (!missing(dens)) {
        if (!is.numeric(dens)) stop("\'dens\' must be numeric")
        if (length(dens) != length(x)) stop("Lengths of \'dens\' and \'x\' must be the same")
        de <-  list(x=sort(x), y=dens[order(x)])
    }
    else { 
        method <- match.arg(method)
        if (method=="frequency") {
            if (missing(nx)) nx <- nclass.Sturges(x)
            xcuts <- seq(min(x), max(x), length=nx)
            dens <- table(cut(x, xcuts)) / length(x)
            de <- list(x=sort(x), y=dens[findInterval(x, xcuts,rightmost.closed=TRUE)][order(x)])
        }
        else if (method=="kernel") de <- density(x, ...)
    }
    if (missing(ny)) ny <- nclass.Sturges(de$y)
    ycuts <- seq(0, max(de$y), length=ny+1)
    rgbmax <- col2rgb(colmax, alpha=TRUE)
    rgbmin <- if (colmin=="transparent") c(col2rgb(colmax, alpha=FALSE), 0) else col2rgb(colmin, alpha=TRUE)
    if (gamma <= 0) stop("gamma must be greater than 0")
    p <- seq(0,1,length=ny-1) ^ {1 / gamma}
    cols <- rgb(rgbmax[1] + p*(rgbmin[1] - rgbmax[1]),
                rgbmax[2] + p*(rgbmin[2] - rgbmax[2]),
                rgbmax[3] + p*(rgbmin[3] - rgbmax[3]),
                alpha=rgbmax[4] + p*(rgbmin[4] - rgbmax[4]),
                maxColorValue=255)
    for (i in 2:ny){
        ## draw rectangles for each region of x with density greather than cut-off
        ## one for each contiguous block in ind 
        ind <- which(de$y >= ycuts[i])
        ints <- seqToIntervals(ind) # from R.utils
        for (j in seq(length=nrow(ints))) {
            if (horiz) 
                rect.fn(xleft=de$x[ints[j,1]],
                     xright=de$x[ints[j,2]],
                     ybottom=at + (i-1)*offset - (!up.left)*width,
                     ytop=at + (i-1)*offset + up.left*width,
                     col=cols[i-1], border = (if(i==ny) cols[i-2] else NA) )
            else
                rect.fn(ybottom=de$x[ints[j,1]],
                     ytop=de$x[ints[j,2]],
                     xleft=at - (i-1)*offset - up.left*width,
                     xright=at - (i-1)*offset + (!up.left)*width,
                     col=cols[i-1], border = (if(i==ny) cols[i-2] else NA) )
        }
    }
    invisible()
}

cistrip <- function(x, at, d, horiz=TRUE, pch = 16, cex=1, lattice=FALSE, ...)
{
    if (missing(d)) d <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 60
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (is.vector(x)) x <- matrix(x, ncol=3)
    n <- nrow(x)
    if (length(at) != n) stop("length of \'at\' should equal the number of estimates in \'x\'")
    pch <- rep(pch, length.out=n)
    cex <- rep(cex, length.out=n)
    if (lattice) {
        points.fn <- panel.points; seg.fn <- panel.segments;
    }
    else {
        points.fn <- points; seg.fn <- segments;
    }
    for (i in 1:n) { 
        if (horiz) { 
            points.fn(x[i,1], at[i], pch=pch[i], cex=cex[i], ...)
            seg.fn(x[i,2], at[i], x[i,3], at[i], ...)
            seg.fn(x[i,2], at[i]-d/2, x[i,2], at[i]+d/2, ...)
            seg.fn(x[i,3], at[i]-d/2, x[i,3], at[i]+d/2, ...)
        }
        else {
            points.fn(at[i], x[i,1], pch=pch[i], cex=cex[i], ...)
            seg.fn(at[i], x[i,2], at[i], x[i,3], ...)
            seg.fn(at[i]-d/2, x[i,2], at[i]+d/2, x[i,2], ...)
            seg.fn(at[i]-d/2, x[i,3], at[i]+d/2, x[i,3], ...)
        }
    }
    invisible()
}

vwstrip <- function(x, dens, at, width, horiz=TRUE, scale=1, limits=c(-Inf,Inf),
                    col="gray", border=NULL, lwd, lty,
                    ticks=NULL, tlen=1, twd, tty, lattice=FALSE, ...) 
{
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(dens)) {
        ## x assumed to be a sample from a distribution, density estimated
        de <- density(x, ...)
        x <- de$x; dens <- de$y
    }
    else {
        if (!is.numeric(dens)) stop("\'dens\' must be numeric")
        if (length(dens) != length(x)) stop("Lengths of \'dens\' and \'x\' must be the same")
        dens <- dens[order(x)]
        x <- sort(x)
    }
    if (lattice) {
        poly.fn <- panel.polygon; seg.fn <- panel.segments
        default.width <- diff(current.panel.limits()[[if(horiz) "ylim" else "xlim"]]) / 20
        default.lwd <- default.twd <- trellis.par.get("add.line")$lwd
        default.lty <- default.tty <- trellis.par.get("add.line")$lty
        default.border <- trellis.par.get("add.line")$col
    }
    else {
        poly.fn <- polygon; seg.fn <- segments
        default.width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 20
        default.lwd <- default.twd <- par("lwd")
        default.lty <- default.tty <- par("lty")
        default.border <- par("fg")
    }
    dens <- dens[x>limits[1] & x<limits[2]]
    x <- x[x>limits[1] & x<limits[2]]
    if (missing(width)) width <- default.width
    if (missing(lwd)) lwd <- default.lwd
    if (missing(lty)) lty <- default.lty
    if (missing(twd)) twd <- default.twd
    if (missing(tty)) tty <- default.tty
    if (is.null(border)) border <- default.border
    dens <- dens / max(dens) * width/2 * scale
    xx <- c(x, rev(x))
    yy <- c(dens, -rev(dens)) + at 
    if (horiz) poly.fn(xx, yy, col=col, border=border, lty=lty, lwd=lwd)
    else poly.fn(yy, xx, col=col, border=border, lty=lty, lwd=lwd)
    if (!is.null(ticks)) {
        dticks <- dens[findInterval(ticks, x)]
        if (horiz) { tx0 <- tx1 <- ticks; ty0 <- at - dticks*tlen; ty1 <- at + dticks*tlen }
        else { ty0 <- ty1 <- ticks; tx0 <- at - dticks*tlen; tx1 <- at + dticks*tlen }
        seg.fn(tx0, ty0, tx1, ty1, lwd=twd, lty=tty)
    }
    invisible()
}

vwstrip.normal <- function(mean, sd, log=FALSE, nx=1000, ...){    
    x <-
        if (log) qlnorm(ppoints(nx, a=0), mean, sd)
        else qnorm(ppoints(nx, a=0), mean, sd)
    dens <- if (log) dlnorm(x, mean, sd) else dnorm(x, mean, sd)
    vwstrip(x=x, dens=dens, ...)
}

bpstrip <- function(x, prob, at, width, horiz=TRUE, scale=1, limits=c(-Inf,Inf),
                    col="gray", border=NULL, lwd, lty,
                    ticks=NULL, tlen=1, twd, tty, lattice=FALSE)
{
    if (!is.numeric(x)) stop("\'x\' must be numeric")
    if (missing(prob)) {
        x <- sort(x)
        prob <- ecdf(x)(x)
    }
    else {
        if (!is.numeric(prob)) stop("\'prob\' must be numeric")
        if (length(prob) != length(x)) stop("Lengths of \'prob\' and \'x\' must be the same")
        prob <- prob[order(x)]
        x <- sort(x)
    }
    if (lattice) {
        poly.fn <- panel.polygon; seg.fn <- panel.segments
        default.width <- diff(current.panel.limits()[[if(horiz) "ylim" else "xlim"]]) / 20
        default.lwd <- default.twd <- trellis.par.get("add.line")$lwd
        default.lty <- default.tty <- trellis.par.get("add.line")$lty
        default.border <- trellis.par.get("add.line")$col
    }
    else {
        poly.fn <- polygon; seg.fn <- segments
        default.width <- diff(par("usr")[if(horiz) 3:4 else 1:2]) / 20
        default.lwd <- default.twd <- par("lwd")
        default.lty <- default.tty <- par("lty")
        default.border <- par("fg")
    }
    prob <- pmin(prob, 1 - prob)
    prob <- prob[x>limits[1] & x<limits[2]]
    x <- x[x>limits[1] & x<limits[2]]
    if (missing(width)) width <- default.width
    if (missing(lwd)) lwd <- default.lwd
    if (missing(lty)) lty <- default.lty
    if (missing(twd)) twd <- default.twd
    if (missing(tty)) tty <- default.tty
    if (is.null(border)) border <- default.border
    prob <- prob / max(prob) * width/2 * scale
    xx <- c(x, rev(x))
    yy <- c(prob, -rev(prob)) + at 
    if (horiz) poly.fn(xx, yy, col=col, border=border, lty=lty, lwd=lwd)
    else poly.fn(yy, xx, col=col, border=border, lty=lty, lwd=lwd)
    if (!is.null(ticks)) {
        pticks <- prob[findInterval(ticks, x)]
        if (horiz) { tx0 <- tx1 <- ticks; ty0 <- at - pticks*tlen; ty1 <- at + pticks*tlen }
        else { ty0 <- ty1 <- ticks; tx0 <- at - pticks*tlen; tx1 <- at + pticks*tlen }
        seg.fn(tx0, ty0, tx1, ty1, lwd=twd, lty=tty)
    }
    invisible()
}

panel.denstrip <- function(...)
{
    denstrip(..., lattice=TRUE)
}

panel.denstrip.normal <- function(...)
{
    denstrip.normal(..., lattice=TRUE)
}

panel.denstrip.legend <- function(...)
{
    denstrip.legend(..., lattice=TRUE)
}

panel.cistrip <- function(...)
{
    cistrip(..., lattice=TRUE)
}

panel.vwstrip <- function(...)
{
    vwstrip(..., lattice=TRUE)
}

panel.vwstrip.normal <- function(...)
{
    vwstrip.normal(..., lattice=TRUE)
}

panel.bpstrip <- function(...)
{
    bpstrip(..., lattice=TRUE)
}

panel.sectioned.density <- function(...)
{
    sectioned.density(..., lattice=TRUE)
}
