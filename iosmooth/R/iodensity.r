iodensity <- function (x, bw, kernel = c("Trap", "Rect", "SupSm"), n.points = 100, x.points)
{
	kernel <- match.arg(kernel, c("Trap", "Rect", "SupSm"))

	if( missing(bw) ) {
		bw = bwadap.numeric(x)
	}
	
	if( missing(x.points) ) {
		range.x <- range(x)
		len <- range.x[2] - range.x[1]
		x.points <- seq(range.x[1]-0.05*len, range.x[2]+0.05*len, len = n.points)
	} else {
        n.points <- length(x.points)
        x.points <- sort(x.points)
    }

	if(kernel=="SupSm") {
		fhat <- ssden(x, bw)
		yp <- fhat(x.points)
   	} else {
   		kernel <- get(paste0("kernel",kernel))
   		yp <- sapply(x.points, function(xpt) mean(kernel( (x-xpt)/bw )/bw))
   	}
   	list(x = x.points, y = yp)
}