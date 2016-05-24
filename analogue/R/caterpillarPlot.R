`caterpillar` <- function(x, ...) {
    UseMethod("caterpillarPlot")
}

`caterpillarPlot` <- function(x, ...) {
    UseMethod("caterpillarPlot")
}

`caterpillarPlot.data.frame` <- function(x, env, useN2 = TRUE, xlab,
                                         ...) {
    ## compute the optima
    opt <- optima(x = x, env = env)
    ## and tolerances
    tol <- tolerance(x = x, env = env, useN2 = useN2)

    if(missing(xlab)) {
        ## grab xlab from env
        xlab <- deparse(substitute(env))
    }

    ## do the plot
    caterpillarPlot.default(x = opt, tol = tol, xlab = xlab, ...)
}

`caterpillarPlot.default` <- function(x, tol, mult = 1, decreasing = TRUE,
                                      labels, xlab = NULL, pch = 21, bg = "white",
                                      col = "black", lcol = col, lwd = 2,
                                      frame.plot = FALSE, ...) {
    ## reorder
    opt <- x
    ord <- order(opt, decreasing = decreasing)
    opt <- opt[ord]
    tol <- tol[ord]

    ## par
    op <- par(yaxs = "i")
    on.exit(par(op))

    ## number of species
    nspp <- length(opt)
    yvals <- seq_len(nspp)

    ## labels == spp names
    if(missing(labels)) {
        labels <- names(opt)
        if(is.null(labels))
            labels <- paste0("Var", yvals)
    }
    linch <- if (!is.null(labels))
        max(strwidth(labels, "inch"), na.rm = TRUE)
    nmai <- par("mai")
    nmai[2L] <- nmai[4L] + linch + 0.1
    par(mai = nmai)

    ## xlab
    if(is.null(xlab))
        xlab <- deparse(substitute(env))

    ## tolerance range
    upr <- opt + (tol * mult)
    lwr <- opt - (tol * mult)

    ## do the plot
    plot(c(lwr, upr), rep.int(yvals, 2), type = "n", axes = FALSE,
         ylab = "", xlab = xlab, ylim = range(0, yvals + 1),
         frame.plot = frame.plot, ...)
    abline(h = yvals, lty = 1, lwd = 0.5, col = "lightgray")
    segments(lwr, yvals, upr, yvals, col = lcol, lwd = lwd, ...)
    points(opt, yvals, pch = pch, bg = bg, col = col, ...)
    axis(side = 1, ...)
    axis(side = 2, labels = labels, at = yvals, las = 1, ...)

    ## return object
    out <- data.frame(Optima = opt, Tolerance = tol)
    invisible(out)
}

`caterpillarPlot.wa` <- function(x, type = c("observed","model"), ...) {
    ## which type of tolerances
    type <- match.arg(type)
    ## extract the optima and tolerances
    opt <- x$wa.optima
    tol <- if(isTRUE(all.equal(type, "observed"))) {
        x$tolerances
    } else {
        x$model.tol
    }
    ## do the plot
    caterpillarPlot.default(x = opt, tol = tol, ...)
}
