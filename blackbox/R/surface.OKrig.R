# adding contours line in extrapolated area
# requires additional member obj$zextrap => modif de predict_surface
# bw is indicator for b&w plot
surface.OKrig <- function (obj, grid.list = NA, extrap = FALSE, graphics.reset = NULL,
                          xlab = NULL, ylab = NULL, main = NULL, zlab = "ln(L)", zlim = NULL,
                          plot.axes={},
                          levels = NULL, type = "C", nx = 80, ny = 80, bw=F, xticks, yticks, ...) {
  out.p <- predict_surface(obj, grid.list = grid.list, extrap = extrap,
                           nx = nx, ny = ny, drop.Z = TRUE)
  if (!is.null(ylab))
    out.p$ylab <- ylab
  if (!is.null(xlab))
    out.p$xlab <- xlab
  if (!is.null(zlab))
    out.p$zlab <- zlab
  if (!is.null(main))
    out.p$main <- main
  xlabori <- xlab;ylabori <- ylab ## ...ori for testing
  rangeindic <- length(unique(which(!is.na(out.p$z))))
  isolinestuff <- function() {
    if (xlabori=="twoNm" & ylabori=="g") {
      isoline(blackbox.getOption("rosglobal")$latt2Ns2)
    }
  }
  if (rangeindic>1) { ## something to plot within the convex hull
    if (type=="p") {
      xlab <- userunit(xlab, format="expression")
      ylab <- userunit(ylab, format="expression")
      ## some of the ticks marks must fall within the plotted range otherwise the wireframe plot fails...
      gnourf <- wireframe(x=out.p$z, row.values=out.p$x, column.values=out.p$y,
                        scales=list(arrows=F, x=xticks, y=yticks, tck=0.5), main=main,
                        drape=T, col.regions=spaMM.colors(64, redshift=1/2), cuts=63,
                        xlab=list(xlab, rot=57), ylab=list(ylab, rot=-25), zlab=zlab,
                        par.settings = list(axis.line = list(col = "transparent")),
                        screen=list(x=-90, y=60, x=30), distance=1/3)
      print(gnourf)
    } else if (type=="C") {
      xlab <- userunit(xlab, format="expression")
      ylab <- userunit(ylab, format="expression")
      if (bw) { ## the number  of levels of the col argument must depend on the local 'levels' (not 'nlevels') variable
        nnul <- 64
        nul <- pretty(range(out.p$z, finite=TRUE), nnul)
        lnul <- length(nul)
        filled.contour(out.p$x, out.p$y, out.p$z,
                       ## il faut aussi caser isoline avant la legende...
                       plot.axes={plot.axes;isolinestuff();
                                  contour(out.p$x, out.p$y, out.p$zextrap, nlevels=10, add=T,
                                          labcex=blackbox.getOption("graphicPars")$labcex)
                         }, ## plot.axes works like a fn definition
                       levels =nul ,
                       col=colorpanel(lnul-1, "grey60", "white"),
                       xlab=xlab, ylab=ylab, main=main, ...)
      } else { ## Color contour plot
        filled.contour(out.p$x, out.p$y, out.p$z,
                       ## il faut aussi caser isoline avant la legende...
                       plot.axes={plot.axes;isolinestuff();
                                  contour(out.p$x, out.p$y, out.p$zextrap, nlevels=10, add=T,
                                          labcex=blackbox.getOption("graphicPars")$labcex)
                         }, ## plot.axes works like a fn definition
                       color.palette=function(v) {spaMM.colors(v, redshift=1/2)}, ## cf p plot above, mais pure fn.
                       nlevels=64, ## gradation of colors; note that the 'contour' call itself has a distinct argument, with default value 10
                       xlab=xlab, ylab=ylab, main=main, ...)
      }
      ## les traits dans la legende sont inherent a filled.contour->rect(). image() (comme dans image.plot() est plus joli)
    }
  } else { ## nothing in convex hull. Trying to produce some ghostly plot...
    if(rangeindic==0) {message.redef(paste("No point in convex envelope for plot (", xlab, ", ", ylab, ")"))}
    if(rangeindic==1) {message.redef(paste("Only one point in convex envelope for plot (", xlab, ", ", ylab, ")"))}
    if (type=="p") {
      xlab <- userunit(xlab, format="expression")
      ylab <- userunit(ylab, format="expression")
      gnourf <- wireframe(x=out.p$z, row.values=out.p$x, column.values=out.p$y,
                        scales=list(arrows=F, x=xticks, y=yticks, tck=0.5), main=main,
                        drape=T, col.regions=rainbow(64, s=0.1), cuts=63, xlab=list(xlab, rot=57), ylab=list(ylab, rot=-25), zlab=zlab,
                        par.settings = list(axis.line = list(col = "transparent")),
                        screen=list(x=-90, y=60, x=30), distance=1/3)
      print(gnourf)
    } else if (type=="C") {
      xlab <- userunit(xlab, format="expression") ## FR two lines corrected (?) from unicode, 040713
      ylab <- userunit(ylab, format="expression")
      filled.contour(out.p$x, out.p$y, out.p$zextrap,
                     ## il faut aussi caser isoline avant la legende...
                     plot.axes={plot.axes;isolinestuff();
                                contour(out.p$x, out.p$y, out.p$zextrap, nlevels=10, add=T,
                                        labcex=blackbox.getOption("graphicPars")$labcex)
                       }, ## plot.axes works like a fn definition
                     color.palette=function(n) {rainbow(n, s=0.1)}, ## ghostly pastel
                     nlevels=64, ## gradation of colors; note that the 'contour' call itself has a distinct argument, with default value 10
                     xlab=xlab, ylab=ylab, main=main, ...)
    }
  }
  invisible()
}
