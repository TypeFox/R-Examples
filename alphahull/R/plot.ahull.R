plot.ahull <-
function (x, add = FALSE, do.shape = FALSE, wlines = c("none", 
    "both", "del", "vor"), wpoints = TRUE, number = FALSE, col = NULL, 
    xlim = NULL, ylim = NULL, lwd = NULL, ...) 
{
    wlines <- match.arg(wlines)
    if (is.null(class(x)) || class(x) != "ahull") {
        cat("Argument is not of class ahull.\n")
        return(invisible())
    }
    if (is.null(col)) {
        col <- c(1, 1, 1, 1, 1, 1)
    }
    else {
        col <- rep(col, length.out = 6)
    }
    if (is.null(lwd)) {
        lwd <- c(1, 1, 2)
    }
    else {
        lwd <- rep(lwd, length.out = 3)
    }
    wlines <- match.arg(wlines)
    plot.dd <- switch(wlines, none = TRUE, both = FALSE, del = FALSE, 
        vor = FALSE)
    if (do.shape) {
        plot.ashape(x$ashape.obj, add = add, wlines = wlines, 
            wpoints = wpoints, number = number, col = col[2:6], 
            xlim = xlim, ylim = ylim, lwd = lwd[1:2], ...)
    }
    else {
        if (plot.dd) {
            if (!add) {
                if (is.null(xlim)) 
                  xlim <- range(x$ashape.obj$x[, 1])
                if (is.null(ylim)) 
                  ylim <- range(x$ashape.obj$x[, 2])
                plot(0, 0, type = "n", xlim = xlim, ylim = ylim, 
                  axes = FALSE, ...)
                axis(side = 1)
                axis(side = 2)
            }
            if (wpoints) {
                points(x$ashape.obj$x, col = col[3], ...)
            }
            if (number) {
                xoff <- 0.02 * diff(range(x$ashape.obj$x[, 1]))
                yoff <- 0.02 * diff(range(x$ashape.obj$x[, 2]))
                text(x$ashape.obj$x[, 1] + xoff, x$ashape.obj$x[, 
                  2] + yoff, 1:(dim(x$ashape.obj$x)[1]), col = col[6], 
                  ...)
            }
        }
        else {
            plot.delvor(x$ashape.obj$delvor.obj, add = add, wlines = wlines, 
                wpoints = wpoints, number = number, col = col[3:6], 
                lwd = lwd[1], xlim = xlim, ylim = ylim, ...)
        }
    }
    arcs <- which(x$arcs[, 3] > 0)
    if (length(arcs) > 0) {
        for (i in arcs) {
            arc(x$arcs[i, 1:2], x$arcs[i, 3], x$arcs[i, 4:5], 
                x$arcs[i, 6], col = col[1], lwd = lwd[3])
        }
    }
    points <- which(x$arcs[, 3] == 0)
    if (length(points) > 0) {
        for (i in points) {
            points(x$arcs[i, 1], x$arcs[i, 2], col = col[1], 
                pch = 19)
        }
    }
}
