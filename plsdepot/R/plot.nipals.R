#'@title Plot NIPALS basic results
#'
#'@description
#'Plot method for objects of class \code{"nipals"}. This function plots either
#'the variables or the observations, on the selected components (i.e. scores). 
#'Variables are plotted inside the circle of correlations. Observations are plotted
#'on a scatter plot.
#'
#'@details
#'Variables are displayed using the correlations in \code{$cor.xt}.
#'
#'@param x An object of class \code{"nipals"}.
#'@param what What to plot. Options are \code{"variables"} and \code{"observations"} 
#'@param comps An integer vector of length two to indicate which components to plot
#'@param cex Character expansion for labels and points.
#'@param col.labels Color for labels of variables.
#'@param pos Position for the labels text (see \code{\link{par}}).
#'@param offset When \code{pos} is specified, this value gives the offset of the labels.
#'@param col.arrows Color for the arrows when plotting variables.
#'@param lwd The line width of arrows.
#'@param length Length of the edges of the arrow head (in inches).
#'@param angle Angle from the shaft of the arrow to the edge of the arrow head.
#'@param col.points Color for the points when \code{what="observations"}.
#'@param pch Plotting character symbol to use (see \code{\link{par}}).
#'@param pt.bg Background (fill) color for the points given by \code{pch=21:25}.
#'@param show.names Logical indicating whether to show names of points.
#'Only used when \code{what="observations"}
#'@param xpd Logical for controlling clipping region of labels and names.
#'@param xlab Title for the x axis.
#'@param ylab Title for the y axis.
#'@param main Main title of the plot.
#'@param col.main Color of main title.
#'@param cex.main Character expansion of main title.
#'@param col.axis Color of axis annotations (tick marks and labels).
#'@param show.grid Logical indicating whether to show grid lines.
#'@param col.grid Color of grid lines. Only used when \code{show.grid=TRUE}.
#'@param \dots Further arguments are passed to labels or points.
#'@author Gaston Sanchez
#'@seealso \code{\link{nipals}}
#'@method plot nipals
#'@S3method plot nipals
#'@examples
#'
#'  \dontrun{
#'  # load data climbing ropes
#'  data(ropes)
#'
#'  # apply nipals with 3 components
#'  nip1 = nipals(ropes[,-1], comps=3)
#'
#'  # plot variables (correlations)
#'  plot(nip1)
#'
#'  # plot observations
#'  plot(nip1, what="obs")
#'
#'  # plot observations with names
#'  plot(nip1, what="obs", show.names=TRUE)
#'  }
#'
plot.nipals <- 
function(x, what = "variables", comps = c(1,2), 
         cex = 1, col.labels = "#5592e3", pos=NULL, offset=0.1,
         col.arrows = "#5b9cf255", lwd=3.5, length = 0, angle = 0,
         col.points = "#5592e3", pch = 21, pt.bg = "#5b9cf255",
         show.names = FALSE, xpd = TRUE, xlab = NULL, ylab = NULL,
         main = NULL, col.main = "gray35", cex.main=1.2,
         col.axis = "gray40",
         show.grid = TRUE, col.grid = "gray95", ...)
{
  # =======================================================
  # check inputs
  # =======================================================
  if (!is.numeric(comps) || length(comps)>2)
    comps = 1:2
  # how many components
  nc = ncol(x$scores)
  if (!comps[1] %in% 1:nc || !(comps[2] %in% 1:nc)) 
    stop("\nInvalid vector of components")
  k1 = comps[1]
  k2 = comps[2]
  if (is.null(xlab)) xlab = paste("axis", k1)
  if (is.null(ylab)) ylab = paste("axis", k2)
  # how many variables
  p = nrow(x$cor.xt)

  ## PLOT
  if (what == "variables")
  { 
    # =======================================================
    # circle of correlations
    # =======================================================
    z = seq(0, 2*pi, l=100)
    op = par(oma = rep(0,4), mar = c(4,3,3,2), pty = "s", xpd=xpd)
    plot.new()
    plot.window(xlim=c(-1.1,1.1), ylim=c(-1.1,1.1), asp=1)
    axis(side=1, labels=FALSE, lwd=0, lwd.ticks=1, col="gray75")
    mtext(seq(-1, 1, 0.5), side=1, at=seq(-1, 1, 0.5), line=0.5, 
          col=col.axis, cex=0.8)
    axis(side=2, labels=FALSE, lwd=0, lwd.ticks=1, col="gray75")
    mtext(seq(-1, 1, 0.5), side=2, at=seq(-1, 1, 0.5), line=0.7, 
          col=col.axis, cex=0.8, las=2)    
    box(col="gray70")
    if (show.grid)
      abline(h=seq(-1, 1, 0.5), v=seq(-1, 1, 0.5), col=col.grid, xpd=FALSE)
    # axis labs
    mtext(xlab, side=1, at=0, line=2.2, col=col.axis)
    mtext(ylab, side=2, at=0, line=2.2, col=col.axis)
    lines(cos(z), sin(z), lwd=2, col="gray90")
    segments(-1, 0, 1, 0, col="gray90")
    segments(0, -1, 0, 1, col="gray90")
    # arrows    
    arrows(x0=rep(0,p), y0=rep(0,p), x1=x$cor.xt[,k1], y1=x$cor.xt[,k2],
           length=length, angle=angle, col=col.arrows, lwd=lwd)
    # text
    if (is.null(pos))
    {
      neg <- x$cor.xt[,k1] < 0
      pos = rep(4, p)
      pos[neg] = 2
    }
    text(x$cor.xt[,k1], x$cor.xt[,k2], labels=rownames(x$cor.xt), 
         col=col.labels, cex=cex, pos=pos, offset=offset, ...)
    if (is.null(main))
      main = "Circle of Correlations"
    #title(main, col.main=col.main)
    mtext(main, side=3, at=-1.2, line=1, adj=0, cex=cex.main, col=col.main)
    par(op)
  } else {
    # =======================================================
    # observations
    # =======================================================
    plot(x$scores[,k1], x$scores[,k2], type="n", axes=FALSE,
         xlab=xlab, ylab=ylab, col.lab=col.axis)
    axis(side=1, col.axis=col.axis, col=col.axis, lwd=0.5, cex.axis=0.9)
    axis(side=2, las=2, col.axis=col.axis, col=col.axis, lwd=0.5, cex.axis=0.9)
    abline(h=0, v=0, col="gray90", lwd=2)
    box(col="gray70")
    if (!show.names)
    {
      points(x$scores[,k1], x$scores[,k2], col=col.points, pch=pch, bg=pt.bg,
           lwd=1.5, cex=cex)
    } else {
      text(x$scores[,k1], x$scores[,k2], labels=rownames(x$scores), 
           col=col.labels, cex=cex, pos=pos, offset=offset, xpd=xpd, ...)
    }
    if (is.null(main))
      main = "Map of Observations"
    #title(main, col.main=col.main)
    mtext(main, side=3, at=par("usr")[1], line=1, adj=0, cex=cex.main, col=col.main)
  }
}
