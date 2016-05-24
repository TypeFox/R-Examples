#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-02 22:11:34 +0100 (ven, 02 dic 2005) $
recurr <- function(series, m, d, start.time=start(series), end.time=end(series), ...) {
	xyz <- embedd(window(series, start=start.time, end=end.time), m=m, d=d)
	D <- dist(xyz)
	D <- as.matrix(D) / max(D)
	grid.x <- time(series)[1:nrow(D)]
	filled.contour(grid.x, grid.x, D, 
		color.palette = function(n) gray(0:(n-1)/(n-1)), 
		xlab = "time", ylab = "time", main = "Recurrence plot", ...)
}
