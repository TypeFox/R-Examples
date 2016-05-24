plot.cpf <- function(x, conf.int = FALSE, mark.time = FALSE,
                     mark = 3, col = 1, lty, ci.lty = 3, xlim, ylim = c(0, 1),
                     xlab = "", ylab = "", bty = "n", legend = TRUE,
                     curvlab = NULL, legend.pos = NULL, legend.bty = "n",
                     ...) {
    
    if (!inherits(x, "cpf"))
        stop("'x' must be of class 'cpf'")
    if (missing(xlim))
        xlim <- c(0, max(x$time))
    ns <- length(x$size.strata)
    if (missing(lty)) {
        lty <- 1:ns
    }
    else {
        if (length(lty) < ns) {
            lty <- lty[1]:ns
        }
    }
    if (length(col) < ns) {
        col <- rep(col[1], ns)
    }

    size <- c(0, x$size.strata)
    levels <- x$X[, ]
    levels[is.null(levels)] <- 1
    who <- rep(levels, x$size.strata)
    
    plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab,
         bty = bty, ...)
    for (i in seq_along(levels)) {
        lines(x$time[who == levels[i]], x$cp[who == levels[i]], type = "s",
              col = col[i], lty = lty[i], ...)
        if (conf.int) {
            lines(x$time[who == levels[i]], x$lower[who == levels[i]],
                  type = "s", lty = ci.lty, col = col[i], ...)
            lines(x$time[who == levels[i]], x$upper[who == levels[i]],
                  type = "s", lty = ci.lty, col = col[i], ...)
        }
        if (mark.time) {
            points(x$time[x$n.event[, 2] != 0][who == levels[i]],
                   x$cp[x$n.event[, 2] != 0][who == levels[i]],
                   pch = mark, col = col[i], ...)
        }
    }
    if (legend && !is.null(x$X)) {
        if (is.null(curvlab)) {
            curvlab <- x$strata
        }
        if (is.null(legend.pos)) {
            legend.pos <- c(xlim[1] + 0.1, ylim[2] - 0.1)
        }
        args <- list(...)
        i <- pmatch(names(args), names(formals("legend")))
        do.call("legend", c(list(legend.pos[1], legend.pos[2],
                                 curvlab, col = col, lty = lty,
                                 bty = legend.bty),
                            args[!is.na(i)]))
    }
}
