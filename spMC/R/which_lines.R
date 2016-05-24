which_lines <-
function(coords, direction, tolerance = pi/8) {
  # Segments indicator
  #
  #     coords coordinates matrix
  #  direction vector (or versor) of choosen direction
  #  tolerance angle tolerance (in radians)

  if (!is.matrix(coords)) coords <- as.matrix(coords)
  n <- dim(coords)[1]
  nc <- dim(coords)[2]
  if (length(direction) != nc) stop("wrong length of direction vector")
  storage.mode(coords) <- "double"
  storage.mode(direction) <- "double"

  id <- .C('wl', n = as.integer(n), nc = as.integer(nc), 
            coords = as.double(coords), dire = as.double(direction),
            tolerance = as.double(tolerance), 
            id = as.integer(vector("integer", n)), PACKAGE = "spMC")$id
  return(as.integer(as.factor(id)))
}

