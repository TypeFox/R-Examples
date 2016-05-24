ladderplot.default <-
function(x, scale=FALSE, col=1, pch=19, lty=1, xlim=c(0.5, ncol(x)+0.5), ylim=range(x), vertical = TRUE, ordered=FALSE, ...)
{
    x <- as.data.frame(x)
    if (scale)
        x <- apply(x, 2, function(x) (x - min(x, na.rm = TRUE))/(max(x, 
            na.rm = TRUE) - min(x, na.rm = TRUE)))
    if (NCOL(x) < 2)
        stop("'x' must have at least 2 columns")
    nr <- nrow(x)
    if (length(col) < nr)
        col <- rep(col, nr)[1:nr]
    if (length(pch) < nr)
        pch <- rep(pch, nr)[1:nr]
    if (length(lty) < nr)
        lty <- rep(lty, nr)[1:nr]
    if (ordered)
        x <- x[,order(colnames(x))]
    y <- data.frame(values=array(unlist(x)), 
        ind=factor(rep(1:ncol(x), each=nrow(x)), labels=colnames(x)))

    id <- match(colnames(x),levels(y$ind))
    if (vertical) {
        with(y, stripchart(values ~ ind, pch=pch, ylim=ylim, xlim=xlim, vertical=vertical, col="white", ...))
        lapply(1:ncol(x), function(i) points(cbind(rep(i,nr), x[,id[i]]), col=col, pch=pch))
        lapply(1:nr, function(i) lines(cbind(id, as.matrix(x)[i,]), col=col[i], lty=lty[i]))
    } else {
        tmp <- xlim
        xlim <- ylim
        ylim <- tmp
        with(y, stripchart(values ~ ind, pch=pch, ylim=ylim, xlim=xlim, vertical=vertical, col="white", ...))
        lapply(1:ncol(x), function(i) points(cbind(x[,id[i]], rep(i,nr)), col=col, pch=pch))
        lapply(1:nr, function(i) lines(cbind(as.matrix(x)[i,], id), col=col[i], lty=lty[i]))
    }
    invisible(NULL)
}

