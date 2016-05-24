## plot method for Procrustes

plot.procr <- function(x, plot.type = "jointplot", plot.dim = c(1,2), main, xlab, ylab, xlim, ylim, asp = 1, pch = 20, 
                       col.X = "cadetblue", col.Y = "gray", col.Yhat = "coral1", 
                       label.conf = list(label = TRUE, pos = 3, cex = 0.8), arrows = TRUE, length = 0.10, ...) 
{
  
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$col)) label.conf$col <- 1
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  labels <- rownames(x$X)
  
  if (plot.type == "jointplot") {
    if (missing(main)) main <- paste("Procrustes Configuration Plot") else main <- main    
    if (missing(xlim)) xlim <- range(c(x$X[,x1], x$Y[,x1]))*1.1
    if (missing(ylim)) ylim <- range(c(x$X[,y1], x$Y[,y1]))*1.1
    
    plot(x$X[,x1], x$X[,y1], main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col.X, ...)
    points(x$Yhat[,x1], x$Yhat[,y1], col = col.Yhat, pch = pch)
    if (label.conf$label) text(x$X[,x1], x$X[,y1], labels = labels, cex = label.conf$cex, pos = label.conf$pos, col = col.X)
    if (label.conf$label) text(x$Yhat[,x1], x$Yhat[,y1], labels = labels, cex = label.conf$cex, pos = label.conf$pos, col = col.Yhat)  
  }
  if (plot.type == "transplot") {
    if (missing(main)) main <- paste("Procrustes Transformation Plot") else main <- main    
    if (missing(xlim)) xlim <- range(c(x$Yhat[,x1], x$Y[,x1]))*1.1
    if (missing(ylim)) ylim <- range(c(x$Yhat[,y1], x$Y[,y1]))*1.1

    plot(x$Y[,x1], x$Y[,y1], main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col.Y)#, ...)
    points(x$Yhat[,x1], x$Yhat[,y1], col = col.Yhat, pch = pch)
    if (label.conf$label) text(x$Y[,x1], x$Y[,y1], labels = labels, cex = label.conf$cex, pos = label.conf$pos, col = col.Y)
    if (label.conf$label) text(x$Yhat[,x1], x$Yhat[,y1], labels = labels, cex = label.conf$cex, pos = label.conf$pos, col = col.Yhat)  
    
    if (arrows) arrows(x$Y[,1], x$Y[,2], x$Yhat[,1], x$Yhat[,2], length = length, col = col.Y, lwd = 0.8)
  }
}

