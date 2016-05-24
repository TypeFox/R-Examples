#' Parse input for coordinates
#'
#' Extract the coordinate locations from the input object.
#'
#' @param x Input object containing the coordinates in some format.
#' @param verbose Print out info of the coordinates.
#' @export

sg_parse_coordinates <- function(x, verbose = FALSE) {

  dunno <- paste("Unable to parse coordinates from input pattern of type:", is(x)[1])

  if(is(x, "ppp")){ # if we have spatstat object
    coord <- cbind(x$x, x$y)
  }
  else if(is(x, "SpatialPoints")) { # if we have sp object
    coord <- as.matrix(x@coords)
  }
  else if(is.data.frame(x)){ # generic dataframe
    coord <- as.matrix(x)
  }
  else if(is.list(x)) { # custom list object
    if(is.null(x$x)) {
      if(is.null(x$coord) | !is.matrix(x$coord)) stop(dunno)
      else {
        coord <- as.matrix(x$coord)
      }
    }
    else if(is.matrix(x$x)){
      coord <- x$x
    }
    else{
      if(length(x$x) != length(x$y)) stop(dunno)
      coord <- cbind(x$x, x$y)
      if(!is.null(x$z)) coord <- cbind(coord, x$z)
    }
  }
  else if(is.matrix(x)){ # just a regular matrix
    coord <- x
  }
  else stop(dunno)
  coord
}
