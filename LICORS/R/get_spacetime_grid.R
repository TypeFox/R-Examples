#' @title 
#' Get an iterator over the space-time coordinates of the field.
#' @description 
#' This function returns a matrix of space-time coordinates of the field.
#' Both for the whole field as well as the truncated field (without the margin)
#' @param dim dimension of the original field (first dimension is time; rest is 
#' space)
#' @param LC.coordinates template of the LC coordinates
#' @export
#' @seealso \code{\link{compute_LC_coordinates}}, \code{\link{setup_LC_geometry}}
#' @examples
#' AA = matrix(rnorm(200), ncol = 10)
#' LC.geom = setup_LC_geometry(speed=1, horizon=list(PLC = 3, FLC = 0), shape = "cone")
#' bb = get_spacetime_grid(dim(AA), LC.geom$coordinates)

get_spacetime_grid <- function(dim, LC.coordinates) {
  
  horizon <- LC_coordinates2control_settings(LC.coordinates)$horizon
  TT <- dim[1]
  space.dim <- as.list(dim[-1])
  if (length(space.dim) > 2) {
    stop("Space dimension > 2 is not implemented.")
  }
  
  truncated.dim <- compute_margin_coordinates(dim, LC.coordinates)
  # initialize with time grid
  grid.spacetime <- expand.grid(time = seq_len(TT))
  grid.spacetime.truncated <- 
    expand.grid(time = horizon$PLC + seq_len(TT - sum(unlist(horizon))))

  if (length(space.dim) == 0) {
    # done
  } else if (length(space.dim) >= 1) {  # add first spatial dimension
    grid.spacetime <- 
      expand.grid(time = grid.spacetime[, "time"], x1 = seq_len(space.dim[[1]]))
    grid.spacetime.truncated <- 
      expand.grid(time = grid.spacetime.truncated[, "time"], 
                  x1 = truncated.dim$space[[1]][1]:truncated.dim$space[[1]][2])
  } else if (length(space.dim) == 2) {
    grid.spacetime <- 
      expand.grid(time = grid.spacetime[, "time"], 
                  x1 = grid.spacetime[, "x1"],
                  x2 = seq_len(space.dim[[2]]))
    grid.spacetime.truncated <- 
      expand.grid(time = grid.spacetime.truncated[, "time"], 
                  x1 = grid.spacetime.truncated[, "x1"],
                  x2 = truncated.dim$space[[2]][1]:truncated.dim$space[[2]][2])
  }
  out <- list(all = as.matrix(grid.spacetime),
              dim.all = dim,
              truncated = as.matrix(grid.spacetime.truncated),
              dim.truncated = truncated.dim$dim
              )
  out$length.all <- prod(out$dim.all)
  out$length.truncated <- prod(out$dim.truncated)
  
  if (length(space.dim) == 0) {
    names(out$dim.all) <- names(out$dim.truncated) <- "time"
  } else {
    names(out$dim.all) <- names(out$dim.truncated) <- 
      c("time", paste0("x", seq_along(space.dim)))
  }
  return(out)
}