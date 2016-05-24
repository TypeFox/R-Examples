#' @title Braces
#' @description Adds curly braces to a plot.
#' 
#' @param xfrom,xto,yfrom,yto start and end points of braces. Direction of brace determined by from and to arguments.
#' @param radius radius of curve in brace.
#' @param col,lty,lwd color, line type, and line width of braces. See \code{\link{par}} for more
#'   details.
#' 
#' @note Orientation of brace is either horizontal or vertical, with axis along largest range of 
#'   x or y in plotting units.
#' 
#' @author Tim Gerrodette \email{tim.gerrodette@@noaa.gov}
#' 
#' @examples
#' plot(x = c(0, 1), y = c(0, 1000), type = "n", xlab= "", ylab = "")
#' braces(xfrom = 0.2, xto = 0.8, yfrom = c(400, 600), yto = c(300, 700))
# 
#' plot(x = c(0, 100), y = c(0, 17), type = "n", xlab = "x", ylab = "y")
#' text(10, 16, "radius =")
#' for (i in 1:8) {
#'   braces(xfrom = 10 * i + 10, xto = 10 * i + 18, yfrom = 1, 
#'          yto = 15, radius = i / 4, lwd = 2)
#'   text(10 * i + 12, 16, round(i / 4, 2))
#' }
# 
#' plot(c(0, 100), c(0, 17), type = "n", xlab = "x", ylab = "y")
#' braces(30, 80, 13, 11, 1)
#' 
#' plot(c(0, 100), c(0, 17), type = "n", xlab = "x", ylab = "y")
#' braces(c(20, 80, 30), c(10,75,40), 1, 15, radius = c(0.2, 0.5, 0.1), 
#'        lwd = c(1, 2, 3), col = 1:2, lty = 1)
#' 
#' plot(c(0, 100), c(0, 17), type = "n")
#' braces(20, 80, 7, 5, 1)
#' braces(20, 80, 13, 15, 1)
#' 
#' @importFrom graphics par segments lines
#' @export
#' 
braces <- function(xfrom, xto, yfrom, yto, radius = 1, col = par("fg"), lty = par("lty"), lwd = par("lwd")) {
 n <- max(length(xfrom), length(xto), length(yfrom), length(yto))
 if (length(xfrom) < n) xfrom <- rep(xfrom, n)[1:n]
 if (length(xto) < n) xto <- rep(xto, n)[1:n]
 if (length(yfrom) < n) yfrom <- rep(yfrom, n)[1:n]
 if (length(yto) < n) yto <- rep(yto, n)[1:n]
 if (length(radius) < n) radius <- rep(radius, n)[1:n]
 if (length(lwd) < n) lwd <- rep(lwd, n)[1:n]
 if (length(lty) < n) lty <- rep(lty, n)[1:n]
 if (length(col) < n) col <- rep(col, n)[1:n]
 
 theta <- seq(0, pi / 2, length = 100)
 sapply(1:length(xfrom), function(i) {
  xmid <- (xfrom[i] + xto[i]) / 2
  ymid <- (yfrom[i] +yto[i]) / 2
  sx <- sign(xto[i] - xfrom[i])
  sy <- sign(yto[i] - yfrom[i])
  vertical <- abs(xto[i] - xfrom[i]) * par("pin")[1] / (par("usr")[2] - par("usr")[1]) < abs(yto[i] - yfrom[i]) * par("pin")[2] / (par("usr")[4] - par("usr")[3])
  if (vertical) {
    rx <- abs(xfrom[i] - xmid)
    ry <- abs(yfrom[i] - yto[i]) / 10 * radius[i]
    if (min(yfrom[i], yto[i]) + ry > ymid - ry) warning("in braces, radius is too large", call. = FALSE)
    segments(xmid, c(yfrom[i] + sy * ry, ymid + sy * ry), xmid, c(ymid - sy * ry, yto[i] - sy * ry), 
             lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xfrom[i] + sx * rx * sin(theta), yfrom[i] + sy * ry - sy * ry * cos(theta), 
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xfrom[i] + sx * rx * sin(theta), yto[i] - sy * ry + sy * ry * cos(theta),
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xto[i] - sx * rx * sin(theta), ymid - ry + ry * cos(theta), 
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xto[i] - sx * rx * sin(theta), ymid + ry - ry * cos(theta), 
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
  } else {
    rx <- abs(xfrom[i] - xto[i]) / 10 * radius[i]
    ry <- abs(yfrom[i] - ymid)
    if (min(xfrom[i], xto[i]) + rx > xmid - rx) warning("in braces, radius is too large", call. = FALSE)
    segments(c(xfrom[i] + sx * rx, xmid + sx * rx), ymid, c(xmid - sx * rx, xto[i] - sx * rx), ymid, 
             lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xfrom[i] + sx * rx - sx * rx *cos(theta), yfrom[i] + sy * ry * sin(theta),
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xto[i] - sx * rx + sx * rx * cos(theta), yfrom[i] + sy * ry * sin(theta),
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xmid + rx - rx * cos(theta), yto[i] - sy * ry * sin(theta), 
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
    lines(xmid - rx + rx * cos(theta), yto[i] - sy * ry * sin(theta),
          lwd = lwd[i], lty = lty[i], col = col[i]
    )
  }
 })
 invisible(NULL)
}