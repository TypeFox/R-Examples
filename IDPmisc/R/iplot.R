### iplot.R

iplot <- function(x,
                  y = NULL,
                  pixs = 1,
                  zmax = NULL,
                  ztransf = function(x){x},
                  colramp = IDPcolorRamp,
                  cex = par("cex"),
                  main = NULL,
                  d.main = 1,
                  cex.main=par("cex.main"),
                  xlab = NULL,
                  ylab = NULL,
                  cex.lab = 1,
                  legend = TRUE,
                  d.legend = 1,
                  cex.axis = par("cex.axis"),
                  nlab.xaxis = 5,
                  nlab.yaxis = 5,
                  minL.axis = 3,
                  las = 1,
                  border = FALSE,
                  oma = c(5,4,1,0)+0.1,
                  mgp = c(2,0.5,0)*cex.axis,
                  tcl = -0.3,
                  ...)
  ## Produces an image scatter plot of a large 2d dataset.

  ## Authors: Andreas Ruckstuhl, Rene Locher
  ## Version 09-04-08
{
  no.xlab <- is.null(xlab)
  no.ylab <- is.null(ylab)
  no.y <- is.null(y)

  if(no.xlab) xlab <- deparse(substitute(x))
  if(no.ylab) ylab <- deparse(substitute(y))

  if(is.data.frame(x)|is.matrix(x)){
    if(ncol(x)>1) {
      if(!no.y) stop("'y' must be NULL when x is a matrix or data.frame with more than 1 column!\n")
      if (no.xlab) xlab <- colnames(x)[1]
      if (no.ylab) ylab <- colnames(x)[2]
      no.y <- FALSE
      y <- x[,2]
      x <- x[,1]
    } else if(ncol(x)==1) {
      if (no.xlab) xlab <- colnames(x)[1]
      x <- x[,1]
    } else stop("Matrix has no columns!\n")
  }
  if(no.y) {
    ylab <- xlab
    xlab <- "Index"
    y <- x
    x <- 1:length(y)
  }

  xy <- NaRV.omit(data.frame(x,y=y))
  x <- xy$x
  y <- xy$y

  mar <- rep(0,4)
  par(oma=oma, mar=mar)

  x.old <- x
  y.old <- y

  xfac <- is.factor(x.old)
  yfac <- is.factor(y.old)

  if(xfac) {
    x <- as.integer(x.old)
  }
  if(yfac) {
    y <- as.integer(y.old)
  }

  if(!(is.vector(x)&is.vector(y)))
    stop("x must be a vector, matrix or data.frame and y must be a vector if present\n")

  w <- par("cin")[1] * 2.54
  h <- par("cin")[2] * 2.54
  mar.legend <- c(0,0,0,3)

  ## 40% of space for color bar and 60% of it for axis labels
  w.legend <- lcm(7*cex.axis*cex*w)
  h.main <- lcm(cex.main*cex*h)
  d.main <- lcm(d.main*cex.main*cex*h)
  d.legend <- lcm(d.legend*cex.main*cex*h)

  if (!is.null(main) & legend) { ## plot title and legend
      lom <- rbind(c(1,rep(0,2)),
                   rep(0,3),
                   c(3,0,2))

      lo <- layout(lom,
                   widths=c(1, d.legend, w.legend),
                   heights=c(h.main, d.main, 1),
                   respect=TRUE)

      iplotMain(main, cex.main, cex=cex)
      iplotLegend(colramp=colramp,ncol=zmax,cex.axis=cex.axis,
                  border=border, mar=c(mar[1],0,mar[3],4)*cex.axis,
                  las=las, tcl=tcl, cex=cex)
  } ## plot title and legend

  if (is.null(main) & legend) { ## plot legend only
      lom <- matrix(c(2,0,1),ncol=3)

      lo <- layout(lom,
                   widths=c(1, d.legend, w.legend),
                   heights=1,
                   respect=TRUE)

      iplotLegend(colramp=colramp,ncol=zmax,cex.axis=cex.axis,
                  border=border, mar=c(mar[1],0,mar[3],4)*cex.axis,
                  las=las, tcl=tcl, cex=cex)
  }## plot legend only

  if (!is.null(main) & !legend) { ## plot title only
      lom <- matrix(c(1,0,2),ncol=1)

      lo <- layout(lom,
                   widths=1,
                   heights=c(h.main, d.main, 1),
                   respect=TRUE)
      iplotMain(main, cex.main, cex=cex)
  } ## plot title only


  ## And now data are plotted
  par(cex=cex, cex.axis=cex.axis, las=las, mar=mar*cex.axis,
       mgp=mgp, tcl=tcl, ...)

  ## drawing labels at
  at.x <- pretty(x,n=nlab.xaxis)
  at.y <- pretty(y,n=nlab.yaxis)

  plot(if (xfac) range(at.x)+0.5*c(-1,+1) else range(x),
       if (yfac) range(at.y)+0.5*c(-1,+1) else range(y),
       xlab=xlab,
       ylab=ylab,
       type="n",
       axes=FALSE,
       main=NULL,
       cex.main=cex.main,
       xpd = NA,
       cex.lab = cex.lab)

  if(is.factor(x.old)) {
    xmin <- min(x,na.rm=TRUE)
    xmax <- max(x,na.rm=TRUE)
    at <- seq(xmin, xmax, by=max(floor((xmax-xmin)/(max(nlab.xaxis-1,1))),1))
    axis(1, at=at,
         labels=abbreviate(levels(x.old)[at],minlength=minL.axis))
  } else {
    axis(1, at = at.x)
  }

  if(is.factor(y.old)) {
    ymin <- min(y,na.rm=TRUE)
    ymax <- max(y,na.rm=TRUE)
    at <- seq(ymin, ymax, by=max(floor((ymax-ymin)/(max(nlab.yaxis-1,1))),1))
    axis(2, at=at,
         labels=abbreviate(levels(y.old)[at],minlength=minL.axis))
  } else {
    axis(2, at = at.y)
  }

  zzmax <- Image(x,y,pixs=pixs,zmax=zmax,ztransf=ztransf,
                 colramp=colramp, factors=c(xfac,yfac))
  if(is.null(zmax)) zmax <- zzmax
  zmax <- max(zmax,2)
  box()

  invisible(zmax)
} ## iplot

