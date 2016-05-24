plot.divchain <- function (x, add = FALSE, ...) 
{
    dotargs <- list(...)
    bty <- dotargs$bty
    bxc <- dotargs$boxcol
    dotargs$boxcol <- NULL
    if (!add) {
        rw <- attr(x, "rw")
        plot(0, 0, type = "n", ann = FALSE, axes = FALSE, xlim = rw[1:2], 
            ylim = rw[3:4])
        if(is.null(bty)) bty <- "n"
        box(bty = bty, col=bxc)
        do.call(title,dotargs)
    }
    lapply(1:nrow(x), function(i, x) {
        do.call(segments, c(as.list(unname(x[i, 1:4])), dotargs))
    }, x = x)
    invisible()
}
