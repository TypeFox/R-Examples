plot.mdsbi <- function(x, plot.dim = c(1,2), sphere = TRUE, col = 1, label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8), 
                       vec.conf = list(col = 1, cex = 0.8, length = 0.1), 
                       identify = FALSE, type = "p", pch = 20, 
                       asp = 1, main, xlab, ylab, xlim, ylim, ...) {
  
  ## argument checks
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  
  if (is.null(label.conf$label)) label.conf$label <- TRUE
  if (is.null(label.conf$pos)) label.conf$pos <- 3
  if (is.null(label.conf$col)) label.conf$col <- 1
  if (is.null(label.conf$cex)) label.conf$cex <- 0.8
  if (identify) label.conf$label <- FALSE
  
  if (is.null(vec.conf$length)) vec.conf$length <- 0.1
  if (is.null(vec.conf$col)) vec.conf$col <- 1
  if (is.null(vec.conf$cex)) vec.conf$cex <- 0.8
  
  if (missing(main)) main <- paste("Configuration Plot") else main <- main
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  
  circle <- function(x, y, r, ...) {
    ang <- seq(0, 2*pi, length = 100)
    xx <- x + r * cos(ang)
    yy <- y + r * sin(ang)
    polygon(xx, yy, ...)
  }
  
  x$conf <- x$model$X           ## configurations
  temp <- rbind(x$conf, t(x$coef))
  if (missing(xlim)) xlim <- range(temp[,x1])*1.1
  if (missing(ylim)) ylim <- range(temp[,y1])*1.1
  
  ## configuration plot
  plot(x$conf[,x1], x$conf[,y1], main = main, type = type, xlab = xlab, ylab = ylab, 
       xlim = xlim, ylim = ylim, pch = pch, asp = asp, col = col, ...)
  if (label.conf$label) text(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), 
                             cex = label.conf$cex, pos = label.conf$pos, 
                             col = label.conf$col)
  if ((any(class(x) == "smacofSP")) && (sphere)) {
    radius <- sqrt(x$conf[2,1]^2 + x$conf[2,2]^2)
    circle(0, 0, radius, lty = 2, border = "lightgray")
  }
  if (identify) {
    identify(x$conf[,x1], x$conf[,y1], labels = rownames(x$conf), cex = label.conf$cex, pos = label.conf$pos, 
             col = label.conf$col)
  }
  
  ## add vectors
  abline(h = 0, v = 0, lty = 2, col = "darkgray")
  n <- ncol(x$coef)
  xycoor <- t(x$coef[c(x1, y1), ])
  posvec <- apply(xycoor, 1, sign)[2,] + 2     
  for (i in 1:n) {
    arrows(0, 0, x$coef[x1, i], x$coef[y1, i], length = vec.conf$length, col = vec.conf$col)
    text(x$coef[x1, i], x$coef[y1, i], labels = colnames(x$coef)[i], col = vec.conf$col, pos = posvec[i], cex = vec.conf$cex)
  }
}
