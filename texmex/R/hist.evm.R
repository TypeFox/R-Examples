hist.evmOpt <-
function(x, xlab, ylab, main, ...){
    # Want parameters as a matrix with one row for passing
    # through to family$rng etc.
    a <- t(x$coefficients)
    u <- x$threshold
    if (!is.finite(u)){ u <- min(x$data$y) }

    # FOLLOWING if BLOCK COMMENTED OUT TO ACCOUNT FOR DIFFERENCE
    # BETWEEN GEV AND GPD. MIGHT HAVE TO DO SOMETHING MORE
    # SENSIBLE LATER.
#    if(a[2] < 0){
#         UpperEndPoint <- u - a[1]/a[2]
#    }
#    else {
        UpperEndPoint <- Inf
#    }

    dat <- x$data$y
    dfun <- x$family$density

    h <- hist(dat, plot = FALSE)
    xx <- seq(u, min(UpperEndPoint, max(h$breaks)), length = 100)
    y <- dfun(xx, a, x)

    breaks <- seq(from=min(dat),to=max(dat),len=nclass.Sturges(dat)+1)

    res <- list(dat=dat, dens=cbind(x=xx, y=y), breaks=breaks)
    oldClass(res) <- "hist.evmOpt"
    res
}

plot.hist.evmOpt <- function(x, xlab=NULL, ylab=NULL, main=NULL, ...){

    if (missing(xlab) || is.null(xlab)) xlab <- "Data"
    if (missing(ylab) || is.null(ylab)) ylab <- ""
    if (missing(main) || is.null(main)) main <- "Histogram and density"

    hist(x$dat, prob = TRUE, ylim = c(0, max(x$dens[, 2])),
         xlab=xlab, ylab=ylab, main=main, breaks = x$breaks, ...)
    lines(x$dens[, 1], x$dens[, 2], col = 4)
    rug(x$dat)
    invisible()
}

print.hist.evmOpt <- plot.hist.evmOpt

