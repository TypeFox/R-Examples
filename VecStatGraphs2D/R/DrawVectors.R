DrawVectors <- function (data, SVGf=0)
{
  	plot.new()

	xini = t(data[, 5])
	yini = t(data[, 6])
	xfin = t(data[, 7])
	yfin = t(data[, 8])
	
    xmin = min(xini, xfin)
	xmax = max(xini, xfin)
	ymin = min(yini, yfin)
	ymax = max(yini, yfin)
	
	plot(xini, yini, main = "Vectors", xlim = c(xmin - 5000, xmax), ylim = c(ymin, ymax))
	
	for (i in 1:(length(xini))) {
	    arrows(xini[i], yini[i], xfin[i], yfin[i])
    }
      
}
