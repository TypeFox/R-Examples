### poster.plot.R


poster.plot <- function(x, y=NULL,  type="p", col=col.fg,
                        col.axis=col.fg, col.lab=col.fg, col.fg="blue",
                        col.bg="lavender",col.box="cornsilk",
                        xlim=NULL, ylim=NULL,
                        xlab=NULL, ylab=NULL, main="", cex=1.2,
                        axes=TRUE, ...){
  ## Convenient xyplot with colored background.
  ##
  ## Authors: Andreas Ruckstuhl, refined by Rene Locher
  ## Version: 2005-10-17

  no.xlab <- (is.null(xlab))
  no.ylab <- (is.null(ylab))

  if(no.xlab) xlab <- deparse(substitute(x))
  if(no.ylab) ylab <- deparse(substitute(y))

  xdim <- dim(x)
  ## coercing 1d matrix into vector
  if(!is.null(xdim)&&xdim[2]==1) x <- as.vector(x)
  if(is.null(y)&is.vector(x)) {
    x.old <- x
    x <- 1:length(x)
    y <- x.old
    if(no.ylab) ylab <- xlab
    xlab <- "Index"
  } else if(is.null(y)&length(dim(x))>=2) {
    t.xlab <- colnames(x)[1]
    t.ylab <- colnames(x)[2]
    y<-x[,2]
    x <- x[,1]
    if(no.xlab) xlab <- t.xlab
    if(no.ylab) ylab <- t.ylab
  }

  opar <- par(bg=col.bg, fg=col.fg, col.axis=col.axis, col.lab=col.lab,
              font=2,cex=cex,...)
   on.exit(opar)

  if(is.null(xlim)) xlim <- range(x)
  if(is.null(ylim)) ylim <- range(y)
  plot(xlim, ylim, type="n", axes=FALSE, ann=FALSE)
  x.usr <- par("usr")
  rect(x.usr[1], x.usr[3], x.usr[2], x.usr[4], col=col.box)
  points(x,y, type=type, col=col); box()
  if(axes) {axis(1); axis(2, las=1)}
  title(main=main,xlab=xlab, ylab=ylab)
  invisible()
} ## poster.plot

