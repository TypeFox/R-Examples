#' @include osm-descriptors.R
{}



merge_ways_nodes <- function(ways, nodes) {
  colnames(ways) <- sprintf("w%s", colnames(ways))
  colnames(nodes) <- sprintf("n%s", colnames(nodes))

  m <- match(ways$wref, nodes$nid)

  dat <- cbind(ways, nodes[m, ])
  # dat <- na.omit(dat)
  dat <- dat[!is.na(dat$nlat), ]

  dat$nid <- NULL
  colnames(dat) <- substring(colnames(dat), 2)

  dat
}



#' Plot osmar object
#'
#' Simple plotting functions to visualize \code{\link{osmar}}
#' objects. Note that for more complex plots, we suggest to convert
#' the objects into \code{sp} and use their plotting functionality.
#'
#' @param x An \code{\link{osmar}} object
#' @param way_args A list of parameters for plotting ways
#' @param node_args A list of parameters for plotting nodes
#' @param ... Ignored
#'
#' @method plot osmar
#'
#' @S3method plot osmar
plot.osmar <- function(x,
                       way_args = list(col = gray(0.7)),
                       node_args = list(pch = 19, cex = 0.1, col = gray(0.3)), ...) {
  n <- dim(x)

  if ( n["nodes"] == 0 )
    stop("Object without nodes.")

  add <- FALSE

  if ( n["ways"] > 0 ) {
    do.call(plot_ways, c(list(x), way_args))
    add <- TRUE
  }

  do.call(plot_nodes, c(list(x, add = add), node_args))
}



#' @param add New plot device or plot on existing onde
#' @rdname plot.osmar
#'
#' @export
plot_nodes <- function(x, add = FALSE, ...) {
  stopifnot(is_osmar(x))

  coords <- x$nodes[[1]][, c("lon", "lat")]

  if ( add )
    points(coords, ...)
  else
    plot(coords, ...)
}



#' @param xlab A x-axis label
#' @param ylab A y-axis label
#' @rdname plot.osmar
#'
#' @export
plot_ways <- function(x, add = FALSE, xlab = "lon", ylab = "lat", ...) {
  stopifnot(is_osmar(x))

  dat <- merge_ways_nodes(x$ways[[3]], x$nodes[[1]])

  rlat <- range(dat$lat, na.rm = TRUE)
  rlon <- range(dat$lon, na.rm = TRUE)

  dat <- split(dat, dat$id)

  if ( !add ) {
    plot(1, type = "n", xlim = rlon, ylim = rlat,
         xlab = xlab, ylab = ylab)
  }

  for ( coord in dat ) {
    if ( nrow(coord) >= 2 ) {
      lines(coord[, c("lon", "lat")], ...)
    }
  }
}

