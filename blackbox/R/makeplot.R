makeplot <- function(xvar, yvar, type="C", main, grid.list, plotobject=blackbox.getOption("fitobject"), ...) {
  ## par(xaxt=..., yxaxt=...) does not work like <plot>(..., xaxt..., yaxt....). Use the latter
  miscOptions <- blackbox.getOption("miscOptions")
  tmp <- maketicks(xvar, axis=1, maxticks=blackbox.getOption("graphicPars")$xmaxticks)
  xlabs <- tmp$labels;xat <- tmp$at;xph <- tmp$phantomat
  tmp <- maketicks(yvar, axis=2, maxticks=blackbox.getOption("graphicPars")$ymaxticks)
  ylabs <- tmp$labels;yat <- tmp$at;yph <- tmp$phantomat
  loccex.axis <- blackbox.getOption("graphicPars")$cex.axis
  if (type=="C") { ## contours: xaxt, yaxt are not arguments for perspective plots
    rosglobal <- blackbox.getOption("rosglobal")
    success <- safeSurface.OKrig(plotobject, extrap=("extrapolateOutOfHull" %innc% miscOptions), type=type, xlab=xvar,
                              xaxt="n", yaxt="n", ## suppresses default ticks on axes
                              ## in log scale, and contout plot generally, ca ne marche pas d'essayer de mettre les axes apres l'appel du plot
                              plot.axes={axis(1, at=xat, labels=mantissExp(xlabs), cex.axis=loccex.axis);axis(2, at=yat, labels=mantissExp(ylabs), cex.axis=loccex.axis);
                                         if (!is.null(xph)) axis(1, at=xph, labels=rep("", length(xph)), tcl=-0.1);
                                         if (!is.null(yph)) axis(2, at=yph, labels=rep("", length(yph)), tcl=-0.1);
                                         points(rosglobal$par[[xvar]], rosglobal$par[[yvar]], pch="+", cex=2)},
                              ylab=yvar, main=main, grid.list=grid.list, ...)
  } else { ## perspective
    xticks <- list(at=xat, labels=mantissExp(xlabs)) ## also $distance (undocumented) to control separation from axis...
    yticks <- list(at=yat, labels=mantissExp(ylabs))
    success <- safeSurface.OKrig(plotobject, extrap=("extrapolateOutOfHull" %innc% miscOptions), type=type, xlab=xvar,
                              ylab=yvar, main=main, grid.list=grid.list, xticks=xticks, yticks=yticks, cex.axis=loccex.axis, ...)
  }
  #   if(success) {
  #      if (type=="C" && xvar=="twoNm" && yvar=="g") {
  #         isoline((canonVPoutput)[2])
  #      }
  #   }
}
