#' @title Computes cartesian coordinates and moments of branches and logs 
#'
#' @description A data frame is populated with branch and log masses, along with \eqn{x}, \eqn{y} cartesian coordinates and \eqn{x}, \eqn{y}, and \eqn{z} moments.
#' \eqn{z} coordinates and moments are calculated only if branches height from the ground (and tilt) have been measured in the field.
#'
#' @param object an object of \code{treeData} class
#' @return an object of class \code{vectors} 
#' @seealso \code{\link{getCoordinatesAndMoment}}
#' @export
treeVectors <- function(object) {
  ## vectors data frame is populated
  vectors <- object$fieldData[, colnames(object$fieldData) %in% c("azimuth", "tipD", "biomass", "height", "tilt", "toBePruned")]

  ## computes cartesian coordinates of branch tip and its x, y and z moments are added to vectors slot
  vectors <- mdply(vectors, getCoordinatesAndMoment, branchesCM = object$branchesCM)
  colnames(vectors) <- c("Azimuth", "Distance", "Height", "Tilt", "toBePruned", "Biomass", "x", "y", "mx", "my", "mz")
  class(vectors) <- c("vector", class(vectors))
  return(vectors)
}

#' @title Plots branches and logs
#'
#' @description Plots branches and logs
#'
#' The 2d plot represents branches and logs as vectors pointing inwards.
#' Branches to be pruned are not shown on graph.
#'
#' @param x      vectors object
#' @param y      unused
#' @param txtcol Colour of text labels, defaults to "grey80"
#' @param ...    Arguments to be passed to plot.default
#' @return \code{NULL}
#' @method plot vector
#' @export
#' @importFrom graphics par
#' @importFrom graphics abline
#' @importFrom graphics text
#' @importFrom graphics plot.default
plot.vector <- function(x, y = NULL, txtcol = "grey80", ...) {
  treeVectors <- x[!x$toBePruned,]

  ## plots branch masses
  ## size of points is proportional to branch biomass
  maxPointSize <- 10
  pointSize <- maxPointSize * as.numeric(as.character(treeVectors$Biomass)) / max(treeVectors$Biomass)

  plot.default(treeVectors[c("x", "y")], cex = pointSize, pch = 13, ...)
  abline(h = 0, v = 0, col = "gray70")
  
  # print vector labels
  chw <- par()$cxy[1] 
  text(treeVectors[c("x", "y")] - chw, labels = row.names(treeVectors), adj = 0, cex = 0.8, col = txtcol) 
}
