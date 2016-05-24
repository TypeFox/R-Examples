`plot.core` <-
function(x, y=NULL, type="h",
         main="", sub="",
         xlim=NULL, ylim=NULL,
         xlab = "Core Size (Number of Residues)",
         ylab = "Total Ellipsoid Volume (Angstrom^3)", 
         axes=TRUE, ann=par("ann"),
         col=par("col"),
         ...) {


  if(is.list(x)) {
    len <- x$length ## hack!! fix later
    x <- x$volume
  } else{ len <- rev(1:length(x)) }
  
  xy <- xy.coords(x, y)
  if (is.null(xlim))
    xlim <- range(xy$x[is.finite(xy$x)])
  if (is.null(ylim))
    ylim <- range(xy$y[is.finite(xy$y)])
##  opar <- par(no.readonly=TRUE)
##  on.exit(par(opar))

  plot.new()
  plot.window(xlim, ylim, ...)
  points(xy$x, xy$y, col=col, type=type, ...)

  if (axes) {
    ax.ind <- c(1,seq(10,length(x),by=10))
    axis(1, at=ax.ind, labels = len[ax.ind])
    axis(2)
    box()
  }
  if (ann) {
    if(is.null(xlab))  xlab=xy$xlab
    if(is.null(ylab))  ylab=xy$ylab
    title(main=main, sub=sub, 
          xlab=xlab, ylab=ylab, ...)
  }

}

