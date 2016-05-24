#' Plot a spatial graph
#'
#' Rudimentary plotting.
#'
#' @param x an 'sg' graph object
#' @param data The point pattern object, same as for computing the 'g'
#' @param which Indices of which out-edges to plot. Default: all
#' @param add Add to existing plot? (default: FALSE)
#' @param addPoints Add points? Will be added if add=FALSE
#' @param points.pch point styling
#' @param points.col point styling
#' @param points.cex point styling
#' @param max.edges limit of edges to try to plot, gets very slow at high count. default 1e4
#' @param ... passed to 'lines' function
#'
#' @importFrom graphics abline axis lines par plot points
#' @importFrom grDevices rgb
#' @export

plot.sg <- function(x, data, which=NULL, add=FALSE,
                    addPoints = FALSE, points.pch=1, points.col=1, points.cex=1,
                    max.edges = 1e4,
                    ...) {
  data <- sg_parse_coordinates(data)

  if(is.null(which)) which <- 1:nrow(data)

  if(ncol(data) == 2) {
    if(!add) {
      plot(NA, NA, xlim=range(data[,1]), ylim=range(data[,2]), asp=1, xlab="x", ylab="y")
      addPoints <- TRUE
    }

    # gather edges, could be big
    which <- sort(which)
    e <- x$edges[which]
    nl <- sapply(e, length)
    ab <- cbind(rep(which, times=nl[which]), unlist(e))

    # unique edges
    ok <- !duplicated(t(apply(ab, 1, sort)))
    ab <- ab[ok,]
    if(sum(ok) > max.edges)
      stop(paste0("Trying to plot too many edges (", sum(ok),"), increase max.edges to override."))

    #
    by_i <- split(data.frame(ab), ab[,1])
    sapply(by_i, function(ab) {
      x0<-  data[ab[1,1],1]
      y0<-  data[ab[1,1],2]
      xo <- data[ab[,2],1]
      yo <- data[ab[,2],2]
      x1 <- as.vector( rbind(x0, xo ))
      y1 <- as.vector( rbind(y0, yo ))
      lines(x1, y1, ...)
    })

    if(addPoints)
      points(data[,1], data[,2], pch=points.pch, col=points.col, cex=points.cex)
  }
  #
  if(ncol(data) == 3) null <- plot3.sg(x, data, which, ...)
  if(ncol(data)>3) stop("Plot only for 2 or 23D.")

}

#' Plot 3d graph
#' @param x sg object
#' @param data coordinates
#' @param which points of which out-edges will be plotted
#' @param ... passed to rgl.lines
#'
#' @importFrom rgl rgl.lines
#' @export
plot3.sg <- function(x, data, which, ...) {

  A <- sg2adj(x)$matrix

  n <- ncol(A)

  which <- sort(which)
  e <- x$edges[which]
  nl <- sapply(e, length)
  if(sum(nl)>1e4) stop("can't handle > 10 000 edges.")

  ab <- cbind(rep(which, times=nl[which]), unlist(e))

  # unique edges
  ok <- !duplicated(t(apply(ab, 1, sort)))
  ab <- ab[ok,]
  #
  by_i <- split(data.frame(ab), ab[,1])
  sapply(by_i, function(ab) {
    x0<-  data[ab[1,1],1]
    y0<-  data[ab[1,1],2]
    z0 <- data[ab[1,1],3]
    xo <- data[ab[,2],1]
    yo <- data[ab[,2],2]
    zo <- data[ab[,2],3]
    x1 <- as.vector( rbind(x0, xo ))
    y1 <- as.vector( rbind(y0, yo ))
    z1 <- as.vector( rbind(z0, zo ))
    rgl.lines(x1, y1, z1, ... )
  })
}
