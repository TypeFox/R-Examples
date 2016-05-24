plot.ashape <-
function (x, add = FALSE, wlines = c("none", "both", "del", "vor"), 
    wpoints = TRUE, number = FALSE, col = NULL, xlim = NULL, 
    ylim = NULL, lwd = NULL, ...) 
{
    wlines <- match.arg(wlines)
    if (is.null(class(x)) || class(x) != "ashape") {
        cat("Argument is not of class ashape.\n")
        return(invisible())
    }
    if (is.null(col)) {
        col <- c(1, 1, 1, 1, 1)
    }
    else {
        col <- rep(col, length.out = 5)
    }
    if (is.null(lwd)) {
        lwd <- 1:2
    }
    else {
        lwd <- rep(lwd, length.out = 2)
    }
    wlines <- match.arg(wlines)
    plot.dd <- switch(wlines, none = TRUE, both = FALSE, del = FALSE, 
        vor = FALSE)
    if (plot.dd) {
        if (!add) {
            if (is.null(xlim)) 
                xlim <- range(x$x[, 1])
            if (is.null(ylim)) 
                ylim <- range(x$x[, 2])
            plot(0, 0, type = "n", xlim = xlim, ylim = ylim, 
                axes = FALSE, ...)
            axis(side = 1)
            axis(side = 2)
        }
        if (wpoints) {
            points(x$x, col = col[2], ...)
        }
        if (number) {
            xoff <- 0.02 * diff(range(x$x[, 1]))
            yoff <- 0.02 * diff(range(x$x[, 2]))
            text(x$x[, 1] + xoff, x$x[, 2] + yoff, 1:(dim(x$x)[1]), 
                col = col[5], ...)
        }
    }
    else {
        plot.delvor(x$delvor.obj, add = add, wlines = wlines, 
            wpoints = wpoints, number = number, col = col[2:5], 
            lwd = lwd[1], xlim = xlim, ylim = ylim, ...)
    }
    ashape <- x$edges
    n2 <- dim(ashape)[1]
    if (n2 >= 1) {
        for (i in 1:n2) {
            segments(ashape[i, "x1"], ashape[i, "y1"], ashape[i, 
                "x2"], ashape[i, "y2"], col = col[1], lwd = lwd[2])
        }
    }
}
