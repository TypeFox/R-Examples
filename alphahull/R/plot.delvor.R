plot.delvor <-
function (x, add = FALSE, wlines = c("both", "del", "vor"), wpoints = TRUE, 
    number = FALSE, col = NULL, xlim = NULL, ylim = NULL, ...) 
{
    wlines <- match.arg(wlines)
    if (is.null(class(x)) || class(x) != "delvor") {
        cat("Argument is not of class delvor.\n")
        return(invisible())
    }
    if (is.null(col)) {
        col <- c(1, 1, 1, 1)
    }
    else {
        col <- rep(col, length.out = 4)
    }
    plot.del <- switch(wlines, both = TRUE, del = TRUE, vor = FALSE)
    plot.vor <- switch(wlines, both = TRUE, del = FALSE, vor = TRUE)
    if (!add) {
        if (is.null(xlim)) 
            xlim <- range(x$x[, 1])
        if (is.null(ylim)) 
            ylim <- range(x$x[, 2])
        plot(0, 0, type = "n", xlim = xlim, ylim = ylim, axes = FALSE, 
            ...)
        axis(side = 1)
        axis(side = 2)
    }
    if (wpoints) {
        points(x$x, col = col[1], ...)
    }
    if (number) {
        xoff <- 0.02 * diff(range(x$x[, 1]))
        yoff <- 0.02 * diff(range(x$x[, 2]))
        text(x$x[, 1] + xoff, x$x[, 2] + yoff, 1:(dim(x$x)[1]), 
            col = col[4], ...)
    }
    if (plot.vor) {
        n.edge <- dim(x$mesh)[1]
        for (i in 1:n.edge) {
            lty.vor = ifelse(x$mesh[i, "bp1"] == 1 | x$mesh[i, 
                "bp2"] == 1, 2, 1)
            segments(x$mesh[i, "mx1"], x$mesh[i, "my1"], x$mesh[i, 
                "mx2"], x$mesh[i, "my2"], lty = lty.vor, col = col[3], 
                ...)
        }
    }
    if (plot.del) {
        plot(x$tri.obj, add = TRUE, do.points = FALSE, col = col[2], 
            ...)
    }
}
