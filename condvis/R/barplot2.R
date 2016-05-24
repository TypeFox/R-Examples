# this code is copied from graphics::barplot.default and edited to return some
# info to allow subsequent editing of the plot
barplot2 <- function (height, width = 1, space = NULL, names.arg = NULL, 
    legend.text = NULL, beside = FALSE, horiz = FALSE, density = NULL, 
    angle = 45, col = NULL, border = par("fg"), main = NULL, 
    sub = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    xpd = TRUE, log = "", axes = TRUE, axisnames = TRUE, cex.axis = 
    par("cex.axis"), cex.names = par("cex.axis"), inside = TRUE, plot = TRUE, 
    axis.lty = 0, offset = 0, add = FALSE, args.legend = NULL, 
    ...) 
{
    if (!missing(inside)) 
        .NotYetUsed("inside", error = FALSE)
    if (is.null(space)) 
        space <- if (is.matrix(height) && beside) 
            c(0, 1)
        else 0.2
    space <- space * mean(width)
    if (plot && axisnames && is.null(names.arg)) 
        names.arg <- if (is.matrix(height)) 
            colnames(height)
        else names(height)
    if (is.vector(height) || (is.array(height) && (length(dim(height)) == 
        1))) {
        height <- cbind(height)
        beside <- TRUE
        if (is.null(col)) 
            col <- "grey"
    }
    else if (is.matrix(height)) {
        if (is.null(col)) 
            col <- gray.colors(nrow(height))
    }
    else stop("'height' must be a vector or a matrix")
    if (is.logical(legend.text)) 
        legend.text <- if (legend.text && is.matrix(height)) 
            rownames(height)
    stopifnot(is.character(log))
    logx <- logy <- FALSE
    if (log != "") {
        logx <- length(grep("x", log)) > 0L
        logy <- length(grep("y", log)) > 0L
    }
    if ((logx || logy) && !is.null(density)) 
        stop("Cannot use shading lines in bars when log scale is used")
    NR <- nrow(height)
    NC <- ncol(height)
    if (beside) {
        if (length(space) == 2) 
            space <- rep.int(c(space[2L], rep.int(space[1L], 
                NR - 1)), NC)
        width <- rep_len(width, NR)
    }
    else {
        width <- rep_len(width, NC)
    }
    offset <- rep_len(as.vector(offset), length(width))
    delta <- width/2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    log.dat <- (logx && horiz) || (logy && !horiz)
    if (log.dat) {
        if (min(height + offset, na.rm = TRUE) <= 0) 
            stop("log scale error: at least one 'height + offset' value <= 0")
        if (logx && !is.null(xlim) && min(xlim) <= 0) 
            stop("log scale error: 'xlim' <= 0")
        if (logy && !is.null(ylim) && min(ylim) <= 0) 
            stop("log scale error: 'ylim' <= 0")
        rectbase <- if (logy && !horiz && !is.null(ylim)) 
            ylim[1L]
        else if (logx && horiz && !is.null(xlim)) 
            xlim[1L]
        else 0.9 * min(height, na.rm = TRUE)
    }
    else rectbase <- 0
    if (!beside) 
        height <- rbind(rectbase, apply(height, 2L, cumsum))
    rAdj <- offset + (if (log.dat) 
        0.9 * height
    else -0.01 * height)
    delta <- width/2
    w.r <- cumsum(space + width)
    w.m <- w.r - delta
    w.l <- w.m - delta
    if (horiz) {
        if (is.null(xlim)) 
            xlim <- range(rAdj, height + offset, na.rm = TRUE)
        if (is.null(ylim)) 
            ylim <- c(min(w.l), max(w.r))
    }
    else {
        if (is.null(xlim)) 
            xlim <- c(min(w.l), max(w.r))
        if (is.null(ylim)) 
            ylim <- range(rAdj, height + offset, na.rm = TRUE)
    }
    if (beside) 
        w.m <- matrix(w.m, ncol = NC)
    if (plot) {
        dev.hold()
        opar <- if (horiz) 
            par(xaxs = "i", xpd = xpd)
        else par(yaxs = "i", xpd = xpd)
        on.exit({
            dev.flush()
            par(opar)
        })
        if (!add) {
            plot.new()
            plot.window(xlim, ylim, log = log, ...)
        }
        xyrect <- function(x1, y1, x2, y2, horizontal = TRUE, 
            ...) {
            if (horizontal) 
                rect(x1, y1, x2, y2, ...)
            else rect(y1, x1, y2, x2, ...)
        }
        if (beside) 
            xyrect(rectbase + offset, w.l, c(height) + offset, 
                w.r, horizontal = horiz, angle = angle, density = density, 
                col = col, border = border)
        else {
            for (i in 1L:NC) {
                xyrect(height[1L:NR, i] + offset[i], w.l[i], 
                  height[-1, i] + offset[i], w.r[i], horizontal = horiz, 
                  angle = angle, density = density, col = col, 
                  border = border)
            }
        }
        if (axisnames && !is.null(names.arg)) {
            at.l <- if (length(names.arg) != length(w.m)) {
                if (length(names.arg) == NC) 
                  colMeans(w.m)
                else stop("incorrect number of names")
            }
            else w.m
            axis(if (horiz) 
                2
            else 1, at = at.l, labels = names.arg, lty = axis.lty, 
                cex.axis = cex.names, ...)
        }
        if (!is.null(legend.text)) {
            legend.col <- rep_len(col, length(legend.text))
            if ((horiz & beside) || (!horiz & !beside)) {
                legend.text <- rev(legend.text)
                legend.col <- rev(legend.col)
                density <- rev(density)
                angle <- rev(angle)
            }
            xy <- par("usr")
            if (is.null(args.legend)) {
                legend(xy[2L] - xinch(0.1), xy[4L] - yinch(0.1), 
                  legend = legend.text, angle = angle, density = density, 
                  fill = legend.col, xjust = 1, yjust = 1)
            }
            else {
                args.legend1 <- list(x = xy[2L] - xinch(0.1), 
                  y = xy[4L] - yinch(0.1), legend = legend.text, 
                  angle = angle, density = density, fill = legend.col, 
                  xjust = 1, yjust = 1)
                args.legend1[names(args.legend)] <- args.legend
                do.call("legend", args.legend1)
            }
        }
        title(main = main, sub = sub, xlab = xlab, ylab = ylab, 
            ...)
        if (axes) 
            axis(if (horiz) 
                1
            else 2, cex.axis = cex.axis, ...)
        invisible(list(w.m = w.m, w.l = w.l, w.r = w.r, height = height))
    }
    else list(w.m = w.m, w.l = w.l, w.r = w.r, height = height)
}