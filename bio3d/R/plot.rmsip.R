"plot.rmsip" <-
  function(x, xlab=NULL, ylab=NULL, col=gray(50:0/50), zlim=c(0,1), ...) {
    
    ##opar <- par(no.readonly = TRUE)
    ##on.exit(par(opar))
    
    if(is.null(xlab))
      xlab <- "a"
    if(is.null(ylab))
      ylab <- "b"
    
    image(1:ncol(x$overlap), 1:nrow(x$overlap), x$overlap,
          col=col, zlim=zlim,
          xlab=xlab, ylab=ylab, ...)
    mtext(paste("RMSIP:", round(x$rmsip, 2)), side=3, line=0.5, at=0.5, adj=0, ...)
    
  }
