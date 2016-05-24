# various high-level plots

# last modified 2014-09-04 by J. Fox

Hist <- function(x, groups, scale=c("frequency", "percent", "density"), xlab=deparse(substitute(x)), 
    ylab=scale, main="", breaks="Sturges", ...){
    xlab # evaluate
    scale <- match.arg(scale)
    ylab
    if (!missing(groups)){
        counts <- table(groups)
        if (any(counts == 0)){
            levels <- levels(groups)
            warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
        }
        levels <- levels(groups)
        hists <- lapply(levels, function(level) if (counts[level] != 0)  
            hist(x[groups == level], plot=FALSE, breaks=breaks)
            else list(breaks=NA))
        range.x <- range(unlist(lapply(hists, function(hist) hist$breaks)), na.rm=TRUE)
        n.breaks <- max(sapply(hists, function(hist) length(hist$breaks)))
        breaks. <- seq(range.x[1], range.x[2], length=n.breaks)
        hists <- lapply(levels, function(level) if (counts[level] != 0) 
            hist(x[groups == level], plot=FALSE, breaks=breaks.)
            else list(counts=0, density=0))
        ylim <- if (scale == "frequency"){
            max(sapply(hists, function(hist) max(hist$counts)))
        }
        else if (scale == "density"){
            max(sapply(hists, function(hist) max(hist$density)))
        }
        else {
            max.counts <- sapply(hists, function(hist) max(hist$counts))
            tot.counts <- sapply(hists, function(hist) sum(hist$counts))
            ylims <- tot.counts*(max(max.counts[tot.counts != 0]/tot.counts[tot.counts != 0]))
            names(ylims) <- levels
            ylims
        }
        save.par <- par(mfrow=n2mfrow(sum(counts != 0)), oma = c(0, 0, if (main != "") 1.5 else 0, 0))
        on.exit(par(save.par))
        for (level in levels){
            if (counts[level] == 0) next
            if (scale != "percent") Hist(x[groups == level], scale=scale, xlab=xlab, ylab=ylab, 
                main=paste(deparse(substitute(groups)), "=", level), breaks=breaks., ylim=c(0, ylim), ...)
            else Hist(x[groups == level], scale=scale, xlab=xlab, ylab=ylab, 
                main=paste(deparse(substitute(groups)), "=", level), breaks=breaks., ylim=c(0, ylim[level]), ...)
        }
        if (main != "") mtext(side = 3, outer = TRUE, main, cex = 1.2)
        return(invisible(NULL))
    }
    x <- na.omit(x)
    if (scale == "frequency") hist(x, xlab=xlab, ylab=ylab, main=main, breaks=breaks, ...)
    else if (scale == "density") hist(x, freq=FALSE, xlab=xlab, ylab=ylab, main=main, breaks=breaks, ...)
    else {
        n <- length(x)
        hist(x, axes=FALSE, xlab=xlab, ylab=ylab, main=main, breaks=breaks, ...)
        axis(1)
        max <- ceiling(10*par("usr")[4]/n)
        at <- if (max <= 3) (0:(2*max))/20
        else (0:max)/10
        axis(2, at=at*n, labels=at*100)
    }
    box()
    abline(h=0)
    invisible(NULL)
}

indexplot <- function(x, labels=seq_along(x), id.method="y", type="h", id.n=0, ylab, ...){
    if (missing(ylab)) ylab <- deparse(substitute(x))
    plot(x, type=type, ylab=ylab, xlab="Observation Index", ...)
    if (par("usr")[3] <= 0) abline(h=0, col='gray')
    ids <- showLabels(seq_along(x), x, labels=labels, id.method=id.method, id.n=id.n)
    if (is.null(ids)) return(invisible(NULL)) else return(ids)
}

lineplot <- function(x, ..., legend){
    xlab <- deparse(substitute(x))
    y <- cbind(...)
    m <- ncol(y)
    legend <- if (missing(legend)) m > 1
    if (legend && m > 1) {
        mar <- par("mar")
        top <- 3.5 + m
        old.mar <- par(mar=c(mar[1:2], top, mar[4]))
        on.exit(par(old.mar))
    }
    if (m > 1) matplot(x, y, type="b", lty=1, xlab=xlab, ylab="")
    else plot(x, y, type="b", pch=16, xlab=xlab, ylab=colnames(y))
    if (legend && ncol(y) > 1){
        xpd <- par(xpd=TRUE)
        on.exit(par(xpd), add=TRUE)
        ncols <- length(palette())
        cols <- rep(1:ncols, 1 + m %/% ncols)[1:m]
        usr <- par("usr")
        legend(usr[1], usr[4] + 1.2*top*strheight("x"), 
            legend=colnames(y), col=cols, lty=1, pch=as.character(1:m))
    }
    return(invisible(NULL))
}

plotDistr <- function(x, p, discrete=FALSE, cdf=FALSE, ...){
    if (discrete){
        if (cdf){
            plot(x, p, ..., type="n")
            abline(h=0:1, col="gray")
            lines(x, p, ..., type="s")
        }
        else {
            plot(x, p, ..., type="h")
            points(x, p, pch=16)
            abline(h=0, col="gray")
        }
    }
    else{
        if (cdf){
            plot(x, p, ..., type="n")
            abline(h=0:1, col="gray")
            lines(x, p, ..., type="l")
        }
        else{
            plot(x, p, ..., type="n")
            abline(h=0, col="gray")
            lines(x, p, ..., type="l")
        }
    }
    return(invisible(NULL))
}


plotMeans <- function(response, factor1, factor2, error.bars = c("se", "sd", "conf.int", "none"),
    level=0.95, xlab=deparse(substitute(factor1)), ylab=paste("mean of", deparse(substitute(response))),
    legend.lab=deparse(substitute(factor2)), main="Plot of Means",
    pch=1:n.levs.2, lty=1:n.levs.2, col=palette(), ...){
    if (!is.numeric(response)) stop("Argument response must be numeric.")
    xlab # force evaluation
    ylab
    legend.lab
    error.bars <- match.arg(error.bars)
    if (missing(factor2)){
        if (!is.factor(factor1)) stop("Argument factor1 must be a factor.")
        valid <- complete.cases(factor1, response)
        factor1 <- factor1[valid]
        response <- response[valid]
        means <- tapply(response, factor1, mean)
        sds <- tapply(response, factor1, sd)
        ns <- tapply(response, factor1, length)
        if (error.bars == "se") sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") sds <- qt((1 - level)/2, df=ns - 1, lower.tail=FALSE) * sds/sqrt(ns)
        sds[is.na(sds)] <- 0
        yrange <-  if (error.bars != "none") c( min(means - sds, na.rm=TRUE), max(means + sds, na.rm=TRUE)) else range(means, na.rm=TRUE)
        levs <- levels(factor1)
        n.levs <- length(levs)
        plot(c(1, n.levs), yrange, type="n", xlab=xlab, ylab=ylab, axes=FALSE, main=main, ...)
        points(1:n.levs, means, type="b", pch=16, cex=2)
        box()
        axis(2)
        axis(1, at=1:n.levs, labels=levs)
        if (error.bars != "none") arrows(1:n.levs, means - sds, 1:n.levs, means + sds,
            angle=90, lty=2, code=3, length=0.125)
    }
    else {
        if (!(is.factor(factor1) | is.factor(factor2))) stop("Arguments factor1 and factor2 must be factors.")
        valid <- complete.cases(factor1, factor2, response)
        factor1 <- factor1[valid]
        factor2 <- factor2[valid]
        response <- response[valid]
        means <- tapply(response, list(factor1, factor2), mean)
        sds <- tapply(response, list(factor1, factor2), sd)
        ns <- tapply(response, list(factor1, factor2), length)
        if (error.bars == "se") sds <- sds/sqrt(ns)
        if (error.bars == "conf.int") sds <- qt((1 - level)/2, df=ns - 1, lower.tail=FALSE) * sds/sqrt(ns)
        sds[is.na(sds)] <- 0
        yrange <-  if (error.bars != "none") c( min(means - sds, na.rm=TRUE), max(means + sds, na.rm=TRUE)) else range(means, na.rm=TRUE)
        levs.1 <- levels(factor1)
        levs.2 <- levels(factor2)
        n.levs.1 <- length(levs.1)
        n.levs.2 <- length(levs.2)
        if (length(pch) == 1) pch <- rep(pch, n.levs.2)
        if (length(col) == 1) col <- rep(col, n.levs.2)
        if (length(lty) == 1) lty <- rep(lty, n.levs.2)
        if (n.levs.2 > length(col)) stop(sprintf("Number of groups for factor2, %d, exceeds number of distinct colours, %d."), n.levs.2, length(col))		
        plot(c(1, n.levs.1 * 1.4), yrange, type="n", xlab=xlab, ylab=ylab, axes=FALSE, main=main, ...)
        box()
        axis(2)
        axis(1, at=1:n.levs.1, labels=levs.1)
        for (i in 1:n.levs.2){
            points(1:n.levs.1, means[, i], type="b", pch=pch[i], cex=2, col=col[i], lty=lty[i])
            if (error.bars != "none") arrows(1:n.levs.1, means[, i] - sds[, i],
                1:n.levs.1, means[, i] + sds[, i], angle=90, code=3, col=col[i], lty=lty[i], length=0.125)
        }
        x.posn <- n.levs.1 * 1.1
        y.posn <- sum(c(0.1, 0.9) * par("usr")[c(3,4)])
        text(x.posn, y.posn, legend.lab, adj=c(0, -.5))
        legend(x.posn, y.posn, levs.2, pch=pch, col=col, lty=lty)
    }
    invisible(NULL)
}
