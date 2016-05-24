cv.select <- function(x, y, degree, interval = NULL, ...) 
{ 		
	if (is.null(interval)) {
		rangex   <- diff(range(x))
		meshx    <- rangex / (length(x) - 1)
		interval <- c( ifelse(degree < 2, meshx / 2, meshx), rangex / 2)
	}
	
	optimize(cv.score, interval, x = x, y = y, degree = degree, ...)$minimum 
}
