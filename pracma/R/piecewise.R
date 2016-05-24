##
##  p i e c e w i s e . R  Piecewise Linear Function
##


piecewise <- function(x, y, abs = FALSE)
{
	n <- length(x)
	areas <- 0.0
	zeros <- if (y[1] == 0) c(x[1]) else c()
	for (i in 2:n) {
		if (y[i]*y[i-1] >= 0) {
			if (y[i] == 0) zeros <- c(zeros, x[i])
			areas <- c(areas, (y[i]+y[i-1]) * (x[i]-x[i-1]) / 2.0)
		} else {
			x0 <- (x[i-1]*y[i] - x[i]*y[i-1])/(y[i] - y[i-1])
			zeros <- c(zeros, x0)
			areas <- c(areas, y[i-1]*(x0-x[i-1])/2.0, y[i]*(x[i]-x0)/2.0)
		}
	}
	area <- if (abs) sum(abs(areas)) else sum(areas)
	return(list(area=area, zeros=zeros))
}
