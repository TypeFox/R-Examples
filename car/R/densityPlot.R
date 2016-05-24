# checked in 2013-06-05 by J. Fox
# 2014-09-04: J. Fox: empty groups produce warning rather than error

densityPlot <- function(x, ...){
    UseMethod("densityPlot")
}

densityPlot.default <- function (x, g, bw="SJ", adjust=1,
    kernel = c("gaussian", "epanechnikov", "rectangular", 
               "triangular", "biweight", "cosine", "optcosine"),
    xlab=deparse(substitute(x)), ylab="Density", 
    col=palette(), lty=seq_along(col), lwd=2, grid=TRUE,
    legend.location="topright", legend.title=deparse(substitute(g)), show.bw=FALSE,
    rug=TRUE, ...) {
    ylab
    if (!is.numeric(x)) stop("argument x must be numeric")
    kernel <- match.arg(kernel)
    if (missing(g)) {
        density <- density(x, bw=bw, adjust=adjust, kernel=kernel)
        if (show.bw) xlab <- paste(xlab, " (bandwidth = ", format(density$bw), ")", sep="")
        plot(density, xlab=xlab, ylab=ylab, main="", type="n", ...)
        if (rug) rug(x)
        if (grid) grid()
        lines(density, lwd=lwd)
    }
    else {
        if (!is.factor(g)) stop("argument g must be a factor")
        counts <- table(g)
        if (any(counts == 0)){
            levels <- levels(g)
            warning("the following groups are empty: ", paste(levels[counts == 0], collapse=", "))
            g <- factor(g, levels=levels[counts != 0])
        }
        legend.title
        valid <- complete.cases(x, g)
        x <- x[valid]
        g <- g[valid]
        levels <- levels(g)
        if (length(bw) == 1) bw <- rep(bw, length(levels))
        if (length(adjust) == 1) adjust <- rep(adjust, length(levels))
        if (length(bw) != length(levels)) stop("number of entries in bw be 1 or must equal number of groups")
        if (length(adjust) != length(levels)) stop("number of entries in adjust must be 1 or must equal number of groups")
        densities <- vector(length(levels), mode="list") 
        names(bw) <- names(adjust) <- names(densities) <- levels
        for (group in levels){
            densities[[group]] <- density(x[g == group], bw=bw[group], 
                                          adjust=adjust[group], kernel=kernel)
        }
        range.x <- range(unlist(lapply(densities, function(den) range(den$x))))
        max.y <- max(sapply(densities, function(den) max(den$y)))
        plot(range.x, c(0, max.y), xlab=xlab, ylab=ylab, type="n", ...)
        if (grid) grid()
        for (i in 1:length(levels)){
            lines(densities[[i]]$x, densities[[i]]$y, lty=lty[i], col=col[i], lwd=lwd)
        }
        if (show.bw){
            bws <- sapply(densities, function(den) den$bw)
            legend <- paste(levels, " (bw = ", format(bws), ")", sep="")
        }
        else legend <- levels
        legend(legend.location, legend=legend, col=col[1:length(levels)], 
               lty=lty, title=legend.title, inset=0.02)
        abline(h=0, col="gray")
        if (rug){
            for (i in 1:length(levels)) rug(x[g == levels[i]], col=col[i])
        }
    }
    return(invisible(NULL))
}

densityPlot.formula <- function(formula, data=NULL, subset, na.action=NULL, xlab, ylab, ...){
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m$xlab <- m$ylab <- m$... <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    if (missing(ylab)) ylab <- "Density"
    if (length(formula) == 3){
        response <- attr(attr(mf, "terms"), "response")
        if (missing(xlab)) xlab <- names(mf)[response]
        g <- mf[, -response]
        densityPlot(mf[[response]], g, xlab=xlab, ylab=ylab, legend.title=names(mf)[-response], ...)
    }
    else if (length(formula) == 2){
        if (missing(xlab)) xlab <- names(mf)
        densityPlot(mf[[1]], xlab=xlab, ylab=ylab,  ...)
    }
    else stop("improper densityPlot formula")   
}
