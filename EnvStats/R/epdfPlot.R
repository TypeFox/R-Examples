epdfPlot <-
function (x, discrete = FALSE, density.arg.list = NULL, plot.it = TRUE, 
    add = FALSE, epdf.col = "black", epdf.lwd = 3 * par("cex"), 
    epdf.lty = 1, curve.fill = FALSE, curve.fill.col = "cyan", 
    ..., type = ifelse(discrete, "h", "l"), main = NULL, xlab = NULL, 
    ylab = NULL, xlim = NULL, ylim = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n.x <- length(x)
    x <- sort(x)
    if (discrete) {
        f.x <- as.numeric(table(x)/n.x)
        x <- unique(x)
    }
    else {
        density.list <- do.call("density", args = c(list(x = x), 
            density.arg.list))
        x <- density.list$x
        f.x <- density.list$y
    }
    names(x) <- NULL
    names(f.x) <- NULL
    if (plot.it) {
        if (!add) {
            if (is.null(main)) 
                main <- paste("Empirical PDF of", data.name)
            if (is.null(xlab)) 
                xlab <- data.name
            if (is.null(ylab)) 
                ylab <- "Relative Frequency"
            if (is.null(xlim)) 
                xlim <- range(x)
            if (is.null(ylim)) 
                ylim <- c(0, max(f.x))
            plot(x, f.x, type = "n", ..., xlim = xlim, ylim = ylim, 
                xlab = xlab, ylab = ylab, main = main)
            arg.list <- list(x = x, y = f.x)
            arg.list <- c(arg.list, checkGraphicsPars(...)$gen.gp.list, 
                list(type = type, col = epdf.col, lwd = epdf.lwd, 
                  lty = epdf.lty))
            do.call("lines", arg.list)
        }
        else lines(x, f.x, ..., type = type, col = epdf.col, 
            lwd = epdf.lwd, lty = epdf.lty)
        if ((!discrete) && curve.fill) {
            n <- length(f.x)
            polygon(c(x, rev(x)), c(f.x, rep(0, n)), border = FALSE, 
                col = curve.fill.col)
        }
    }
    invisible(list(x = x, f.x = f.x))
}
