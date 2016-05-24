#' @title Setup light cone geometry
#' @description 
#' \code{setup_LC_geometry} sets up the light cone geometry for LICORS.
#'
#' @param horizon a list with \code{PLC} and \code{FLC} horizon
#' @param speed speed of propagation
#' @param space.dim dimension of the spatial grid. Eg. \code{2} if the data is 
#' a video ( = image sequences).
#' @param shape shape of light cone: \code{'cone'}, \code{'tube'}, 
#' or \code{'revcone'}.
#' @keywords manip
#' @return 
#' A list of class \code{"LC"}.
#' 
#' @export
#' @seealso \code{\link{LC-utils}}, \code{\link{compute_LC_coordinates}}
#' @examples
#' aa = setup_LC_geometry(horizon = list(PLC = 3, FLC = 1), speed = 1, 
#'                        space.dim = 1, shape = "cone")
#' aa
#' plot(aa)
#' summary(aa)
#' 

setup_LC_geometry <- function(horizon = list(PLC = 1, FLC = 0), 
                              speed = 1, space.dim = 1, shape = "cone") {
  
  if (is.null(horizon$FLC)) {
    horizon$FLC <- 0
  }
  
  out <- list(horizon = horizon,
              speed = speed,
              space.dim = space.dim,
              shape = shape)
  
  out$coordinates <- 
    list(PLC = compute_LC_coordinates(horizon = horizon$PLC, speed = speed, 
                                      space.dim = space.dim, shape = shape, 
                                                  type = "PLC"),
         FLC = compute_LC_coordinates(horizon = horizon$FLC, speed = speed, 
                                      space.dim = space.dim, shape = shape, 
                                      type = "FLC"))
  out$n.p <- nrow(out$coordinates$PLC)
  out$n.f <- nrow(out$coordinates$FLC)
  class(out) <- "LC"
  return(out)
} 
