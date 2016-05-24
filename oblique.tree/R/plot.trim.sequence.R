plot.trim.sequence <- function(x, ..., type = "l", ylim = range(x$dev),
    order = c("decreasing", "increasing"))
{
    if(missing(type) && inherits(x, "trim")) type <- "S"
    if(is.null(x$method)) x$method <- "deviance"
    order <-  match.arg(order)
    if(order == "increasing")
        sign <- +1 else sign <- -1
    plot(sign*x$comp, x$dev, axes = FALSE,
         xlab = "comp", ylab = x$method,
         type = type, ylim = ylim, ...)
    box()
    axis(2, ...)
    xaxp <- par("xaxp")
    pos <- sign*seq(xaxp[1], xaxp[2], diff(xaxp[-3])/xaxp[3])
    if(pos[1] == 0) pos[1] <- 1
    n <- length(pos)
    maxsize <- max(x$comp)
    if(pos[n] > maxsize) pos[n] <- maxsize
    axis(1, at = sign*pos, labels = pos, ...)
    axis(3, at = sign * x$comp, labels = format(signif(x$h, 2)), ...)
    invisible()
}

