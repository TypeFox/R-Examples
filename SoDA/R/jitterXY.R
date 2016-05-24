jitterXY <- function(x = c(0,1), y = c(0,1),
               xscale = 1, yscale = 1) {
	xlim <- range(x, na.rm = TRUE, finite = TRUE)
	ylim <- range(y, na.rm = TRUE, finite = TRUE)
	plot.window(xlim, ylim)
	cxy <- par("cxy")
	if(missing(y))
	   x + xscale * cxy[[1]] * runif(length(x))
	else if(missing(x))
	   y + yscale * cxy[[2]] * runif(length(y))
	else
	   list(x = x + xscale * cxy[[1]] * runif(length(x)),
	        y = y + yscale * cxy[[2]] * runif(length(y)))
}

jitterX <- function(x, scale = 1)
   jitterXY(x, xscale = scale)
jitterY <- function(y, scale = 1)
   jitterXY(y = y, yscale = scale)
