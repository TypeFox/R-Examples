#' Plot 3D surface of a template brain
#'
#' @details This function will work immediately for the standard
#'   \link{templatebrain} defined in the package documentation. If passed an
#'   object called e.g. FCWB it expects to find another object named FCWB.surf
#'   containing the surface information. If you follow this naming convention
#'   for user-defined refbrains it will work for them as well.
#'
#' @param x the template brain to plot.
#' @param col the color of the surface.
#' @param alpha the alpha value of the surface.
#' @param ... extra arguments to pass to \code{\link[rgl]{plot3d}}.
#' @export
#' @method plot3d templatebrain
#' @importFrom rgl plot3d
plot3d.templatebrain <- function(x, col='grey', alpha=0.3, ...) {
  plot3d(get(paste0(deparse(substitute(x)), ".surf")), col=col, alpha=alpha, ...)
}
