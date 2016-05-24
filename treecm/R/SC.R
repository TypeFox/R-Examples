#' @title Compute the slenderness ratio
#'
#' @description slenderness ratio is an important index of stability of trees and branches
#'
#' @note The coefficient takes into account branch angle: 
#' \eqn{SL_c=\frac{L}{D} \cdot (1 + cos \alpha)}, 
#' where \eqn{\alpha} is the branch angle (0 degrees = horizontal, 90 degrees vertical),
#' \eqn{L} is branch length in m, \eqn{D} is branch diameter in cm
#' Vertical branches have \eqn{SL = SL_c}
#'
#' @param x the data frame holding the measures needed to perform the computation
#' @param diameter The name of the data frame column holding diameter of the branch
#' @param length The name of the data frame column holding length of the branch
#' @param tilt The name of the data frame column holding tilt of the branch
#' @return slenderness ratio
#' @references Mattheck, C. and Breloer, H. \emph{The Body Language of Trees: A Handbook for Failure Analysis (Research for Amenity Trees)} 1995, HMSO (London)
branchSR <- function(x, diameter, length, tilt) {
  tiltRad <- as.double(x[tilt]) * pi / 180
  SR <- as.double(x[length]) / as.double(x[diameter]) * 100
  round(SR * ( 1 + cos(tiltRad)), digits = 0)
}

#' @title Computes slenderness ratio for tree branches
#'
#' @description slenderness ratio is an important index of stability of trees and branches
#'
#' @note The coefficient takes into account branch angle: 
#' \eqn{SL_c=\frac{L}{D} \cdot (1 + cos(a))}, 
#' where \eqn{a} is the branch angle (0 degrees = horizontal, 90 degrees vertical). 
#' Vertical branches have \eqn{SL = SL_c}
#'
#' @param treeObject an object of \code{treeData} class
#' @param vectorObject an object of \code{vectors} class
#' @return an object of class \code{SR}
#' @references Mattheck, C. and Breloer, H. \emph{The Body Language of Trees: A Handbook for Failure Analysis (Research for Amenity Trees)} 1995, HMSO (London)
#' @export
#' @seealso \code{\link{branchSR}}
treeSR <- function(treeObject, vectorObject) {
  SR <- treeObject$fieldData[, colnames(treeObject$fieldData) %in% c("dBase", "length", "tilt")]
  SR <- cbind(SR, vectorObject$Azimuth)
  SR <- cbind(SR, apply(SR, 1, branchSR, diameter = "dBase", length   = "length", tilt     = "tilt"))
  colnames(SR) <- c("diameter", "length", "tilt", "azimuth", "SR")
  class(SR) <- c("SR", class(SR))
  return(SR)
}


#' @title Plots slenderness ratio of branches
#'
#' @description Plots the branches as arrows whose length is proportional to their slenderness ratio.
#' A red circle holds ``safe'' branches (\eqn{SR_c \leq 70}).
#'
#' @note A circle is drawn to encompass 
#' the 70- values for slenderness ratio. Branches with 70+ values for the slenderness
#' ratio are considered dangerous. Please note that Mattheck coefficient is corrected to account 
#' for branch tilt (the more it deviates from the verticality the higher its coefficient) 
#'
#' @param x      SR object
#' @param y      unused
#' @param safeSR SR threshold, risky branches are red-coloured
#' @param ...    Arguments to be passed to plot.default
#' @return \code{NULL}
#' @method plot SR
#' @seealso \code{\link{treeSR}}
#' @export
#' @importFrom graphics par
#' @importFrom graphics arrows
#' @importFrom graphics text
#' @importFrom graphics lines
#' @importFrom graphics plot
plot.SR <- function(x, y = NULL, safeSR = 70, ...) {
  
  Circle <- function(t, a) {
    a * cos(t) + 1i * a * sin(t)
  }
  
  xy <- toCartesianXY(x$azimuth, x$SR)
  xyL <- length(xy)
  xyCoord <- data.frame(cbind(xy[1:(xyL/2)], xy[((xyL/2)+1):xyL]))
  colnames(xyCoord) <- c("x", "y")
  rownames(xyCoord) <- rownames(x)
  
  plot(xyCoord, type="n", asp = 1, ...)
  chw <- par()$cxy[1] 
  
  xyCoord <- cbind(xyCoord, x$SR)
  colnames(xyCoord) <- c("x", "y", "SR")
  safe   <- xyCoord[xyCoord$SR <= safeSR, colnames(xyCoord) %in% c("x", "y")]
  unsafe <- xyCoord[xyCoord$SR > safeSR, colnames(xyCoord) %in% c("x", "y")]
  
  if (nrow(safe) > 0) {
    arrows(0, 0, safe$x, safe$y, col = "green")
    text(safe$x - chw, safe$y - chw, labels = row.names(safe), adj = 0, cex = 0.8, col = "green") 
  }
  if (nrow(unsafe) > 0) {
    arrows(0, 0, unsafe$x, unsafe$y, col = "red")
    text(unsafe$x - chw, unsafe$y - chw, labels = row.names(unsafe), adj = 0, cex = 0.8, col = "red") 
  }
  
  t <- seq(0, 2 * pi, by=0.01)
  center <- 0 + 0i
  lines(center + Circle(t, safeSR), col = "red", lwd = 2)
}

