plotDens2 <- function(mle, col=gray(seq(.9,.3,len=30)), border=NA,  
         key=TRUE, n.key=10, round.key=2, 
         numbers=FALSE, round.numbers=2, cex.numbers=0.6, 
         xlim=NULL, ylim=NULL, zlim=NULL, breaks=NULL, 
         xlab="", ylab="", main="", sub="")
{
    p <- mle$p
    rects <- mle$rects

    # Check input
    if (!is.vector(p)) stop("invalid argument 'p'")
    if (sum(p) > 1+1e-5 | sum(p) < 1-1e-5) warning("sum(p) is not equal to 1")
    if (!is.matrix(rects) | ncol(rects) != 4) {
        stop("invalid argument 'rects': it must be a matrix with 4 columns")
    }
    if (sum(rects[, 1] > rects[, 2]) + sum(rects[, 3] > rects[, 4]) >= 1) {
        stop("invalid argument 'rects': each row x1,x2,y1,y2 must 
satisfy x1<=x2 and y1<=y2")
    }
    if (sum(rects[,1]==rects[,2]) + sum(rects[,3]==rects[,4]) >=1){
        stop("cannot create density plot since some maximal 
intersectsions have equal x-coordinates or equal y-coordinates")
    }
    if (length(p) != nrow(rects)) stop("length(p) must equal nrow(rects)")
    if (!is.null(breaks)) {
        if (length(breaks) != (length(col) + 1)) {
            stop("length(breaks) must equal length(col)+1")
        }
    }
    if (!is.logical(key)) stop("argument 'key' must be logical")
    if (!is.logical(numbers)) stop("argument 'numbers' must be logical")
  
    dens <- p/((rects[, 2] - rects[, 1]) * (rects[, 4] - rects[, 3]))

    # determine colors for density 
    np <- length(p)
    ncol <- length(col)
    densmin <- min(dens)
    densmax <- max(dens)
    if (ncol==1 | densmin==densmax) denscol <- rep(1,np) else {
        if (is.null(breaks)){
            if (!is.null(zlim)){
                densmin <- max(densmin,zlim[1])
                densmax <- min(densmax,zlim[2])
            }
            a <- (densmax-densmin)/(ncol-1)
            breaks <- seq(from=densmin-a/2,by=a,to=densmax+a/2)
        } 
        denscol <- cut(dens, breaks, labels=FALSE)
    }

    if (!key) {
        plotRects(rects, col = col[denscol], density = -1, border = border, 
            xlim = xlim, ylim = ylim, xlab=xlab, ylab=ylab, main=main, sub=sub)
    } else {
        # Determine x axis and y axis of plot
        if (is.null(xlim)) xrange <- range(rects[, 1:2]) else xrange <- xlim
        if (is.null(ylim)) yrange <- range(rects[, 3:4]) else yrange <- ylim
       
        xmin <- xrange[1] - 0.04 * (xrange[2] - xrange[1])
        xmax <- xrange[2] + 0.08 * (xrange[2] - xrange[1])
        ymin <- yrange[1] - 0.04 * (yrange[2] - yrange[1])
        ymax <- yrange[2] + 0.04 * (yrange[2] - yrange[1])

        # Plot maximal intersections
        plot(c(0, 0), type = "n", xlim = c(xmin, xmax), ylim = c(ymin, 
            ymax), xaxs = "i", yaxs = "i", xlab = xlab, ylab = ylab, 
            main = main, sub = sub)
        plotRects(rects, col = col[denscol], density = -1, border = border, 
            add = TRUE)

        # Add the key
        key.x1 <- xmax - 0.04 * (xmax - xmin)
        key.x2 <- xmax
        if (densmin==densmax){
            rect(key.x1, ymin, key.x2, ymax, col=col[1], border="black")
            axis(4, at=mean(c(ymin,ymax)), labels=round(breaks, round.key))
        } else if (ncol==1){
            rect(key.x1, ymin, key.x2, ymax, col=col, border="black")
            axis(4, at=c(ymin,ymax), labels=round(c(densmin,densmax),round.key))
        } else {
            key.y <- ymin + c(0:ncol) * (ymax - ymin)/ncol
            rect(rep(key.x1, ncol), key.y[1:ncol], rep(key.x2, ncol), 
                 key.y[2:(ncol + 1)], col = col, border = NA)
            rect(key.x1, ymin, key.x2, ymax, border = "black")
            if (ncol <= 1.5 * n.key) {
                ind.label <- c(1:length(breaks))
            } else {
                a <- ceiling(log2(ncol/(1.5 * n.key)))
                ind.label <- seq(from = 2^a, by = 2^a, to = ncol+1)
            }
            axis(4, at = key.y[ind.label], labels = round(breaks, 
                 round.key)[ind.label])
        }
    }

    if (numbers) {
        midpoints.x <- (rects[, 1] + rects[, 2])/2
        midpoints.y <- (rects[, 3] + rects[, 4])/2
        text(midpoints.x, midpoints.y, round(p, round.numbers), 
            cex = cex.numbers)
    }
}

