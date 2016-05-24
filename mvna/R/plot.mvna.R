plot.mvna <- function(x, tr.choice, xlab = "Time",
                      ylab = "Cumulative Hazard", col = 1, lty, xlim, ylim,
                      conf.int = FALSE, level = 0.95,
                      var.type = c("aalen", "greenwood"),
                      ci.fun = c("log", "linear", "arcsin"),
                      ci.col = col, ci.lty = 3,
                      legend = TRUE, legend.pos, curvlab, legend.bty = "n",
                      ...) {

    if (!inherits(x, "mvna")) {
        stop("'x' must be of class 'mvna'")
    }
    ref <- paste(x$trans[, 1], x$trans[, 2])
    if (missing(tr.choice)) {
        tr.choice <- ref
    }
    if (sum(tr.choice %in% ref) != length(tr.choice)) {
        stop("Names of the possible transitions and 'tr.choice' must match")
    }

    object <- mvna::summary.mvna(x, level = level, var.type = var.type,
                                 ci.fun = ci.fun)[tr.choice]
    lt <- length(object)

    if (missing(lty)) {
        lty <- seq_len(lt)
    }
    else if (length(lty) < lt) {
        lty <- lty * rep(1, lt)
    }
    if (length(col) < lt)
        col <- col * rep(1, lt)

    if (missing(xlim)) {
        xlim <- c(0, max(sapply(object, function(x) max(x$time))))
    }
    if (missing(ylim)) {
        if (conf.int) {
            ylim <- c(0, max(sapply(object, function(x) max(x$upper))))
        } else {
            ylim <- c(0, max(sapply(object, function(x) max(x$na))))
        }
    }

    ## time to do the plotting
    plot(xlim, ylim, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, type = "n", ...)

    for (i in seq_len(lt)) {
        lines(object[[i]]$time, object[[i]]$na, type = "s",
              col = col[i], lty = lty[i], ...)
    }

    if (conf.int) {
        if (length(ci.col) < lt)
            ci.col <- ci.col * rep(1, lt)
        if (length(ci.lty) < lt)
            ci.lty <- ci.lty * rep(1, lt)
        for (i in seq_len(lt)) {
            lines(object[[i]]$time, object[[i]]$lower, type = "s",
                  col = ci.col[i], lty = ci.lty[i], ...)
            lines(object[[i]]$time, object[[i]]$upper, type = "s",
                  col = ci.col[i], lty = ci.lty[i], ...)
        }
    }

    if (legend) {
        if (missing(legend.pos))
            legend.pos <- "topleft"
        if (missing(curvlab))
            curvlab <- tr.choice
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
        do.call("legend", c(list(xx, yy, curvlab, col=col, lty=lty, bty = legend.bty),
                            args[!is.na(ii)]))
    }
    invisible()
}
