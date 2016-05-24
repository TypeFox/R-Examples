##' Calculate maximum distance from 'home' for each trip
##'
##' This function returns a distance from a given 'home' coordinate for each individual trip. 
##' Use the \code{home} argument to provide a single, common 2-element (x,y or lon,lat) coordinate. If \code{home}
##' is \code{NULL} (the default), then each individual trip's first location is used. 
##' @param x trip object
##' @param home see details
##' @seealso \code{\link[sp]{spDistsN1}}
##' @return numeric vector of distances in km (for longlat), or in the units of the trip's projection 
homedist <- function(x, home = NULL) {
  if (!is.null(home)) {
    if (is.numeric(home) & length(home) == 2) {
    home <- matrix(home, ncol = 2)
    } else {
      stop("stop, home must be a 2-element numeric")
    }
  }  
  ## iterate over individual trips
  tor <- getTORnames(x)
  ids <- unique(x[[tor[2L]]])
  dists <- numeric(length(ids))
  names(dists) <- as.character(ids)
  longlat <- is.projected(x)
  if (is.na(longlat)) {
    longlat <- TRUE
    warning("coordinate system is NA, assuming longlat . . .")
  }
  for (i in seq_along(ids)) {
    x0 <- coordinates(x[x[[tor[2L]]] == ids[i], ])
    if (is.null(home)) home <- x0[1, , drop = FALSE]
    dists[i] <- max(spDistsN1(x0, home))
  }
  dists
}