safeSurface.OKrig <- function(plotobject, grid.list,
                           extrap=("extrapolateOutOfHull" %innc% blackbox.getOption("miscOptions")), ...) {
  surface.OKrig(plotobject, grid.list, extrap=extrap, ...)
  return(invisible(TRUE))
}
