xyplot.tags <-
function (xx, yy, tag, log = NULL, xlim = NULL, ylim = NULL, 
    xlab = deparse(substitute(xx)), ylab = deparse(substitute(yy)), 
    taglab = deparse(substitute(tag)), main = "", ...) 
{
    frame()
    temp.x <- remove.na(cbind(xx, yy))
    x <- temp.x$x[1:temp.x$n, 1]
    y <- temp.x$x[1:temp.x$n, 2]
    if (main == "") 
        if (taglab == "") 
            banner <- ""
        else banner <- paste("Plot of 'values' for", taglab)
    else banner <- main
    tag[is.na(tag)] <- "+"
    if (is.null(log)) 
        log <- ""
    plot(x, y, log = log, xlim = xlim, ylim = ylim, type = "n", 
        xlab = xlab, ylab = ylab, main = banner, ...)
    text(x, y, tag, ...)
    invisible()
}
