#' @title Get LC configuration from a (N+1)D field
#'
#' @description 
#' \code{compute_margin_coordinates} computes the coordinates (boundary) of the 
#' margin of the field.
#' 
#' @param dim a vector with the dimensions of the field (time, space1, space2, ..., spaceN)
#' @param LC.coordinates template of the LC coordinates
#' @keywords manip
#' @export
#' @seealso \code{\link{compute_LC_coordinates}}
#' @examples
#' LC_geom = setup_LC_geometry(speed=1, horizon=list(PLC = 3, FLC = 0), shape = "cone")
#' 
#' data(contCA00)
#' 
#' aa = compute_margin_coordinates(dim(contCA00$observed), LC_geom$coordinates)
#' aa
compute_margin_coordinates <- function(dim, LC.coordinates){
  
  controls <- LC_coordinates2control_settings(LC.coordinates)
  horizon <- c(controls$horizon$PLC, controls$horizon$FLC)
  if (length(horizon) == 1){
    horizon <- c(horizon, 0)
  }
  
  TT <- dim[1]
  space.dim <- as.list(dim[-1])

  out <- list(time = c(horizon[1] + 1, TT - horizon[2]))
  
  out$space <- list()
  out$dim <- c(diff(out$time) + 1)
  for (ss in seq_along(space.dim)) {
    out$space[[ss]] <- c(controls$space.cutoff + 1, 
                         space.dim[[ss]] - controls$space.cutoff )
    out$dim <- c(out$dim, diff(out$space[[ss]]) + 1)
  }
  
  return(out)
}


