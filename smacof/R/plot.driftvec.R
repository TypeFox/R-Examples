plot.driftvec <- function(x, main, xlim, ylim, xlab = "Dimension 1", ylab = "Dimension 2", pch = 20, asp = 1, 
                          col.conf = "black", col.drift = "lightgray", 
                          label.conf = list(label = TRUE, pos = 3, col = "black", cex = 0.8), ...) {
  
  if (missing(main)) main <- "Drift Vectors" else main <- main
  if (missing(xlim)) xlim <- range(c(x$fitsym$conf[,1], x$driftcoor[,1]))*1.1
  if (missing(ylim)) ylim <- range(c(x$fitsym$conf[,2], x$driftcoor[,2]))*1.1
  
  plot(x$fitsym$conf, main = main, asp = asp, col = col.conf, pch = pch, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
  if (label.conf$label) text(x$fitsym$conf[,1], x$fitsym$conf[,2], labels = rownames(x$fitsym$conf), 
                             cex = label.conf$cex, pos = label.conf$pos, col = label.conf$col)
  arrows(x$fitsym$conf[,1], x$fitsym$conf[,2], x$driftcoor[,1], x$driftcoor[,2], length = 0.10, lwd = 0.5, col = col.drift)  
}

