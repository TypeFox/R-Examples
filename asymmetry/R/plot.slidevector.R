plot.slidevector <-
function(x,plot.dim=c(1,2),yplus=0,xlab,ylab, ...) {
  #
  # add defaults
  #
  x1 <- plot.dim[1]
  y1 <- plot.dim[2]
  if (missing(xlab)) xlab <- paste("Dimension", x1,sep = " ") else xlab <- xlab
  if (missing(ylab)) ylab <- paste("Dimension", y1,sep = " ") else ylab <- ylab
  plot(x$confi[,x1], x$confi[,y1], xlab=xlab,ylab=ylab, ...)
  arrows(0,0,x$slvec[1],x$slvec[2])
  if(!is.null(rownames(x$confi))){
   text(x$confi[,x1], yplus+x$confi[,y1], rownames(x$confi))
  } else {
    text(x$confi[,x1], yplus+x$confi[,y1], c(1:x$nobs))
  }
}
