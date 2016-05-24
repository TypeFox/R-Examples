plot.clos.etm <- function(x, xlab = "Time",
                          ylab.e = "Expected LOS", ylab.w = "Weights",
                          xlim, ylim.e, ylim.w, col.e = c(1, 2), col.w = 1,
                          lty.e = c(1, 1), lty.w = 1, legend = TRUE,
                          legend.pos, curvlab, legend.bty = "n", ...) {
    if (!inherits(x, "clos.etm")) {
        stop("'x' must be a 'clos.etm' object")
    }
    if (missing(xlim)) {
        xlim <- c(0, max(x$w.time))
    }
    if (missing(ylim.e)) {
        ylim.e <- c(0, max(c(x$phi.case, x$phi.control)))
    }
    if (missing(ylim.w)) {
        ylim.w <- c(0, max(x$weights))
    }
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    split.screen(figs=matrix(c(rep(0,2), rep(1,2), c(0, 0.6), c(0.7, 1)), ncol=4))
    screen(2)
    op <- par(mar=c(2, 5, 2, 1))
    plot(c(0,x$w.time), c(0, x$weights), type = "s", axes = FALSE, lty = lty.w, xlim = xlim,
         ylim = ylim.w , xlab = xlab , ylab = ylab.w, col=col.w, ...)
    axis(side=2)
    box()
    par(op)
    screen(1)      
    op <- par(mar=c(5, 5, 4, 1))
    plot(x$time, x$phi.case, type = "s", lty = lty.e[1], xlim = xlim,
         ylim = ylim.e, xlab = xlab, ylab = ylab.e, col = col.e[1], ...)
    lines(x$time, x$phi.control, type = "s", lty = lty.e[2], col = col.e[2], ...)
    par(op)
    if (legend == TRUE) {
        if (missing(legend.pos))
            legend.pos <- "bottomright"
        if (missing(curvlab))
            curvlab <- c("Intermediate event by time t", "No intermediate event by time t")
        if (is.list(legend.pos)) legend.pos <- unlist(legend.pos)
        if (length(legend.pos) == 1) {
            xx <- legend.pos
            yy <- NULL
        }
        if (length(legend.pos) == 2) {
            xx <- legend.pos[1]
            yy <- legend.pos[2]
        }
        args <- list(...)
        ii <- pmatch(names(args),
                     names(formals("legend")[-charmatch("bty",names(formals("legend")))]))
        do.call("legend", c(list(xx, yy, curvlab, col = col.e, lty = lty.e, bty = legend.bty),
                            args[!is.na(ii)]))
    }
    close.screen(all.screens = TRUE)
    invisible()
}
