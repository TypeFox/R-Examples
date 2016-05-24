makeplottypes <- function(xvar,
                        yvar,
                        types=c("C", "p"),
                        main, grid.list,
                        x11bool= blackbox.getOption("interactiveGraphics"),
                        ...) {
  if ("p" %in% types) { ## 3-D plot of the likelihood surface
    if (x11bool) provideDevice(bbdefaultPars=TRUE) ## if screen output, one window for _each_ graphic
    makeplot(xvar, yvar, type="p", main=main,
             grid.list= grid.list, ...)
  }
  ## Note that 2nd dimension is fixed to its estimated value, dimension 1 of kriging is assigned to x axis, etc
  ## You can change these settings by following the same syntax (but one value has to be fixed).
  ## By default (x, y) points out of the convex hull of the data are clipped out of the figure.
  ## safeSurface.OKrig(..., extrap=TRUE, ...) will keep them, but then the margins may not be good-looking.
  ## These margins may be then fixed by plotrange(fittedNames[<dim>])[<subset>] (e.g. <subset>=2:79)
  if ("C" %in% types) { ## 2-D contour plot of the likelihood surface
    if (x11bool) provideDevice(bbdefaultPars=TRUE) ## one window for _each_ graphic
    ## 2D Contour plot of the likelihood surface:
    makeplot(xvar, yvar, type="C",
             main=main,
             grid.list= grid.list, ...)
  }
}
