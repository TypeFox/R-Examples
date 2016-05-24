`mefalogo` <-
function(type = c("cont", "fill"), labels = c("str", "r"), adj=c(0,0), scale=c(1,1),
new = TRUE, cex = 1, ann=FALSE, axes=FALSE, xlim=c(0,10)*scale[1]+adj[1], ylim=c(0,10)*scale[2]+adj[2], ...)
{
    if (length(type) == 2) type <- type[1]
    if (length(labels) == 2) labels <- labels[1]
    type <- match.arg(type, c("cont", "fill"))
    labels <- match.arg(labels, c("str", "r"))
    x <- c(5,10,10,5)
    y <- c(0,0,5,5)
    z <- list()
    z[[1]] <- cbind(x, y)
    z[[2]] <- cbind(x-1, y+1)
    z[[3]] <- cbind(x-2, y+2)
    z[[4]] <- cbind(c(0,0,3,3), c(2,7,7,2))
    z[[5]] <- cbind(c(3,8,8,3), c(7,7,10,10))
    for (i in 1:5) {
        z[[i]][,1] <- z[[i]][,1]*scale[1]+adj[1]
        z[[i]][,2] <- z[[i]][,2]*scale[2]+adj[2]}

    if (type == "fill") {
        border <- NA
        col <- c("red", "orange", "yellow", "blue", "green")
    } else {
        border <- NULL
        col <- rep("white", 5)
    }
    if (new) plot(0, type="n", ann=ann, axes=axes, ylim=ylim, xlim=xlim, ...)
    for (i in 1:5) polygon(z[[i]], col=col[i], border=border)
    if (labels == "r") {
        text(5.5*scale[1]+adj[1], 4.5*scale[2]+adj[2], "R", col="blue", cex=15*cex)
    } else {
        text(5.5*scale[1]+adj[1], 4.5*scale[2]+adj[2], "count data matrix\n(x$xtab)", cex=cex)
        text(7.5*scale[1]+adj[1], 0.5*scale[2]+adj[2], "segments\n(x$segm)", cex=cex)
        text(1.5*scale[1]+adj[1], 4.5*scale[2]+adj[2], "data frame\nfor samples\n(x$samp)", cex=cex)
        text(5.5*scale[1]+adj[1], 8.5*scale[2]+adj[2], "data frame for taxa\n(x$taxa)", cex=cex)
    }
invisible()
}


