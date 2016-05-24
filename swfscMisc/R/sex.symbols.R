#' @title Sex Symbols
#' @description Plots male and female symbols on current plot.
#' 
#' @param x,y the \code{x} and \code{y} coordinates on the current plot.
#' @param sex a numeric vector containing the values 1 (male) or 2 (female). If of length one, then
#'   value is recycled for all symbols.
#' @param col,lwd,cex color, line width, and character expansion for each point. \code{lwd} and 
#'   \code{col} are recycled as necessary to cover all points. See \code{\link{par}} for more details.
#'   
#' @author Tim Gerrodette \email{tim.gerrodette@@noaa.gov}
#' 
#' @examples
#' x <- runif(20, 0, 10)
#' y <- runif(20, 0, 200)
#' plot(x, y, type = "n")
#' sex.symbols(x, y, sex = 1:2, cex = 1.5, lwd = c(1.5, 4), col = c("blue", "red"))    
#' 
#' @importFrom graphics par strwidth strheight points arrows segments
#' @export
#' 
sex.symbols <- function(x, y, sex = 1, col = par("fg"), lwd = par("lwd"), cex = 1) {
	if (length(y) != length(x)) warning("In function sex.symbols, length of y positions does not equal length of x positions")
	if (length(sex) != length(x)) sex <- rep(sex, length(x))[1:length(x)]
	if (length(col) != length(x)) col <- rep(col, length(x))[1:length(x)]
	if (length(lwd) != length(x)) lwd <- rep(lwd, length(x))[1:length(x)]
  
  cex <- cex[1]
  if (!all(sex %in% c(1, 2))) warning("sex must be coded 1 or 2")
	xm <- strwidth("m") * cex
	ym <- strheight("m") * cex
  xm.in <- par("cin")[1] * cex
  ym.in <- par("cin")[2] * cex
  
	if (is.vector(x)) {
		for (i in 1:length(x)) {
			if (sex[i] == 1 | sex[i] == 2) points(x[i], y[i], pch = 1, lwd = lwd[i], cex = cex, col = col[i])
			if (sex[i] == 1) arrows(x[i] + 0.2 * xm, y[i] + 0.2 * ym, x[i] + 0.60 * xm, y[i] + 0.60 * ym, lwd = lwd[i], 
                              angle = 30, code = 2, length = max(xm.in, ym.in) / 4, col = col[i])
			if (sex[i] == 2) {
        segments(x[i], y[i] - 0.25 * ym, x[i], y[i] - 0.65 * ym, lwd = lwd[i], col = col[i])
 		    segments(x[i] - 0.25 * xm, y[i] - 0.45 * ym, x[i] + 0.25 * xm, y[i] - 0.45 * ym, lwd = lwd[i], col = col[i])
  		}
  	}
	}
  invisible(NULL)
}