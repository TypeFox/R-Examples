Dotplot <- function(x, by, bin=FALSE, breaks, xlim, xlab=deparse(substitute(x))){
    dotplot <- function(x, by, bin=FALSE, breaks, xlim,
                        xlab=deparse(substitute(x)), main="", correction=1/3, correction.char=1, y.max){
        if (bin) hist <- hist(x, breaks=breaks, plot=FALSE)
        if (missing(by)){
            y <- if (bin) hist$counts else table(x)
            x <- if (bin) hist$mids else sort(unique(x))
            plot(range(x), 0:1, type="n", xlab=xlab, ylab="", main=main, axes=FALSE,
                 xlim=xlim)
            y.limits <- par("usr")[3:4]
            char.height <- correction.char*par("cxy")[2]
            axis(1, pos=0)
            if (missing(y.max)) y.max <- max(y)
            abline(h=0)
            cex <- min(((y.limits[2] - y.limits[1])/char.height)/
                           y.max, 2)
            for (i in 1:length(y)){
                if (y[i] == 0) next
                points(rep(x[i], y[i]), cex*correction*char.height*seq(1, y[i]), pch=16, cex=cex,
                       xpd=TRUE)
            }
            return(invisible(NULL))
        }
        else{
            if (missing(xlim)) xlim <- range(x)
            levels <- levels(by)
            n.groups <- length(levels)
            save.par <- par(mfrow=c(n.groups, 1))
            on.exit(par(save.par))
            if (bin){
                for(level in levels){
                    # compute histograms by level to find maximum count
                    max.count <- 0
                    hist.level <- hist(x[by == level], breaks=hist$breaks, plot=FALSE)
                    max.count <- max(max.count, hist.level$counts)
                }
                for (level in levels){
                    dotplot(x[by == level], xlab=xlab, main=paste(label.by, "=", level),
                            bin=TRUE, breaks=hist$breaks, xlim=xlim, correction=1/2, 
                            correction.char=0.5, y.max=max.count)
                }
            }
            else {
                y <- table(x, by)
                for (level in levels){
                    dotplot(x[by == level], xlab=xlab, main=paste(label.by, "=", level),
                            xlim=xlim, correction=1/2, correction.char=0.5, y.max=max(y))
                }
            }
        }
    }
    if (!is.numeric(x)) stop("x must be a numeric variable")
    if (!missing(by) && !is.factor(by)) stop("by must be a factor")
    force(xlab)
    if (missing(by)){
        x <- na.omit(x)
    }
    else{
        label.by <- deparse(substitute(by))
        keep <- complete.cases(x, by)
        x <- x[keep]
        by <- by[keep]
    }
    if (missing(xlim)) xlim <- range(x)
    force(xlab)
    if (missing(breaks))breaks <- "Sturges"
    if (missing(by)) dotplot(x=x, bin=bin, breaks=breaks, xlim=xlim, xlab=xlab)
    else dotplot(x=x, by=by, bin=bin, breaks=breaks, xlim=xlim, xlab=xlab)
}