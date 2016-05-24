corplot <- function(mat, col=mcmcplotsPalette(11, "sequential"), outline=TRUE, greek = FALSE, legend.scale=0.75, mar=c(5, 4, 1, 1) + 0.1, ...){
    opar <- par(mar=mar)
    on.exit(par(opar))
    labels <- dimnames(mat)[[2]]
    if (greek) {
      labels <- .to.greek(labels)
    }
    p <- nrow(mat)
    mat[lower.tri(mat, diag=TRUE)] <- NA
    mat <- abs(mat[, p:1])
    ## com <- combn(1:p, 2)
    ## x <- com[2]
    ## y <- com[1]
    ## z <- mat[lower.tri(mat)]
    image(mat, axes=FALSE, col=col)
    axis(2, at=0:(p-1)/(p-1), labels=c(rev(labels[-1]), NA), las=2, tick=FALSE, ...)
    axis(1, at=0:(p-1)/(p-1), labels=c(labels[-length(labels)], NA), las=2, tick=FALSE, ...)
    ol.x <- seq(par("usr")[1], par("usr")[2], length.out=p + 1)
    ol.y <- seq(par("usr")[3], par("usr")[4], length.out=p + 1)
    if (outline){
        segments(ol.x[1], ol.y[-c(1, length(ol.x))], rev(ol.x)[-c(1, length(ol.x))], ol.y[-c(1, length(ol.x))], col="black", lty=1)
        segments(ol.x[-c(1, length(ol.x))], rev(ol.y)[-c(1, length(ol.x))], ol.x[-c(1, length(ol.x))], ol.y[1], col="black", lty=1)
    }
    ncat <- length(col)
    ## box.height <- (par("usr")[4] - par("usr")[3]) / (p + 1)
    legend.height <- legend.scale*(par("usr")[4] - par("usr")[3])
    box.height <- legend.height / (ncat + 1)
    box.width <- (par("usr")[2] - par("usr")[1]) / (p + 1)
    xleft <- par("usr")[2] - box.height
    xright <- par("usr")[2]
    ytop <- par("usr")[4] - 0:(ncat-1)*box.height
    ybottom <- par("usr")[4] - 1:ncat*box.height
    ## xleft <- rev(ol.x)[2]
    ## xright <- rev(ol.x)[1]
    ## ytop <- rev(ol.y)[1] - 0:(ncat-1)*box.height
    ## ybottom <- rev(ol.y)[1] - 1:ncat*box.height
    rect(xleft, ybottom, xright, ytop, col=rev(col))
    text(xleft, ybottom + 0.5*box.height, labels=formatC(seq(1, 0, length.out=ncat), digits=1, format="f"), pos=2)
}

