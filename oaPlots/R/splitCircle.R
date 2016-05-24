


#' Function for drawing the bottom right portion of a split circle
#' @param x x location of the circle center
#' @param y y location of the circle center
#' @param radius radius of the circle
#' @param splitAngle angle (in radians) that splits the color in two halves
#' @param nv number of vertices used to draw the circle
#' @param border binary whether to include a border on the circle
#' @param col color for the semi-circle
#' @param lty line type used for drawing the circle polygon
#' @param lwd line width used for darwing the circle polygon
#' @return none, semi-circle is drawn to the current device
#' @author Jason Waddell
#' @noRd
draw.halfCircleBR <- function (x, y, radius, splitAngle, nv = 100, border = NULL, col = NA, lty = 1, 
		lwd = 1) 
{
	angle.inc <- 2 * pi/nv
	angles <- c(seq(0, splitAngle , by = angle.inc), seq(splitAngle+pi, 2*pi, by = angle.inc))
	if (length(col) < length(radius)) 
		col <- rep(col, length.out = length(radius))
	for (circle in 1:length(radius)) {
		xv <- cos(angles) * radius[circle] + x
		yv <- sin(angles) * radius[circle] + y
		polygon(xv, yv, border = border, col = col[circle], lty = lty, 
				lwd = lwd)
	}
	invisible(list(x = xv, y = yv))
}

#' Function for drawing the top left portion of a split circle
#' @param x x location of the circle center
#' @param y y location of the circle center
#' @param radius radius of the circle
#' @param splitAngle angle (in radians) that splits the color in two halves
#' @param nv number of vertices used to draw the circle
#' @param border binary whether to include a border on the circle
#' @param col color for the semi-circle
#' @param lty line type used for drawing the circle polygon
#' @param lwd line width used for darwing the circle polygon
#' @return none, semi-circle is drawn to the current device
#' @author Jason Waddell
#' @noRd
draw.halfCircleTL <- function (x, y, radius, splitAngle, nv = 100, border = NULL, col = NA, lty = 1, 
		lwd = 1) 
{
	angle.inc <- 2 * pi/nv
	angles <- seq(splitAngle, splitAngle+pi, by = angle.inc)
	if (length(col) < length(radius)) 
		col <- rep(col, length.out = length(radius))
	for (circle in 1:length(radius)) {
		xv <- cos(angles) * radius[circle] + x
		yv <- sin(angles) * radius[circle] + y
		polygon(xv, yv, border = border, col = col[circle], lty = lty, 
				lwd = lwd)
	}
	invisible(list(x = xv, y = yv))
}


#' Function for drawing a split circle (two differently colored semicircles)
#' @param x x location of the circle center
#' @param y y location of the circle center
#' @param radius radius of the circle
#' @param splitAngle angle (in radians) that splits the color in two halves
#' @param nv number of vertices used to draw the circle
#' @param border binary whether to include a border on the circle
#' @param col1 color of the first semicircle
#' @param col2 color of the second semicircle
#' @param lty line type used for drawing the circle polygon
#' @param lwd line width used for darwing the circle polygon
#' @return none, split circle is drawn to the current device
#' @author Jason Waddell
#' @export
#' @examples 
#' plot(-1, -1, xlim = c(0, 1), ylim = c(0,1), type = "n")
#' splitCircle(x = 0.5, y = 0.5, radius = 0.48,
#' 		splitAngle = pi/4, nv = 1000, border = NA,
#' 		col1 = "blue", col2 = "red")
splitCircle <- function(x, y, radius, splitAngle = pi/4, nv = 100, border = NA, col1 = NA, col2 = NA,
		lty = 1, lwd = 1){
	
	
	draw.halfCircleTL(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
			col = col1)
	draw.halfCircleBR(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
			col = col2)
	
	segments(x0 = x + cos(splitAngle)*radius, x1 = x - cos(splitAngle)*radius, 
			y0 = y + sin(splitAngle)*radius, y1 = y - sin(splitAngle)*radius, col = "white", lwd = 2)
}




#