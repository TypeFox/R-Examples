`plot.pca.scree` <-
function(x, y=NULL, type="o", pch=18,
         main="", sub="",
         xlim=c(0,20), ylim=NULL, 
         ylab="Proporton of Variance (%)",
         xlab="Eigenvalue Rank",
         axes=TRUE, ann=par("ann"),
         col=par("col"), lab=TRUE,
         ...) {

  if(is.list(x))  x=x$L # output from pca.xyz()
  
  PC <- c(1:length(x))
  percent <- (x/sum(x))*100
  cumv<-cumsum(percent)

  #xy <- xy.coords(x, y)
  xy <- xy.coords(percent, y=NULL)
  if (is.null(xlim))
    xlim <- range(xy$x[is.finite(xy$x)])
  if (is.null(ylim))
    ylim <- range(xy$y[is.finite(xy$y)])
  
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  plot.new()
  plot.window(xlim, ylim)#, ...)
  points(xy$x, xy$y, col=col, type=type, pch=pch)#, ...)


  if (axes) {
    axis(1, at=c(1:8,max(xlim)))
    axis(2, at=c(round(percent[c(1:5)],1),0))
    ##axis(1)
    ##axis(2)
    box()
  }
  if(lab) {
    text(c(1:4), percent[c(1:4)], labels=round(cumv[c(1:4)],1), pos=4)
    text(c(8,20), percent[c(8,20)], labels=round(cumv[c(8,20)],1), pos=3)
  }
  if (ann) {
    if(is.null(xlab))  xlab=xy$xlab
    if(is.null(ylab))  ylab=xy$ylab
    title(main=main, sub=sub, 
          xlab=xlab, ylab=ylab, ...)
  }
  out<-list(pc=PC,percent=percent,cumv=cumv)
}

