#' Display an image slice in 3D
#'
#' @param brain template brain (e.g. IS2) of the slice.
#' @param slice path to the slice PNG image to display.
#' @param ... extra arguments to pass to \code{\link[rgl]{persp3d}}.
#' @export
display_slice <- function(brain, slice, ...) {
  UseMethod("display_slice")
}

#' @method display_slice default
display_slice.default <- function(brain, slice, ...) {
}

#' @importFrom rgl persp3d
#' @method display_slice templatebrain
display_slice.templatebrain <- function(brain, slice, z=NULL, ...) {
  bbox <- brain$BoundingBox

  x <- c(bbox[1, 1], bbox[2, 1])
  y <- c(bbox[1, 2], bbox[2, 2])
  if(is.null(z)) {
    z <- (bbox[2, 3] - bbox[1, 3]) / 2
  }
  z <- matrix(rep(z, 4), nrow=2)

  persp3d(x,y,z, col='white', specular='black', texture=slice, xlab='', ylab='', zlab='', box=F, axes=F, add=T, ...)
}
