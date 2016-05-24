ioksmooth <- function (x, y, bw = 0.10, kernel = c("Trap", "Rect", "SupSm"), n.points = 100, x.points) 
{
	kernel <- match.arg(kernel, c("Trap", "Rect", "SupSm"))
	
    if (missing(y)) 
        stop("y must be supplied.\nFor density estimation use iodensity()")
	if (missing(x.points)) {
		range.x <- range(x)
		len <- range.x[2] - range.x[1]
		x.points <- seq(range.x[1]-0.05*len, range.x[2]+0.05*len, len = n.points)
	} else {
        n.points <- length(x.points)
        x.points <- sort(x.points)
    }    
    
	if(kernel=="SupSm") {
 		numer <- ssden(x, y, bw)
 		denom <- ssden(x, 1, bw)
 		yp <- numer(x.points)/denom(x.points)
    } else {
		kernel <- get(paste0("kernel",kernel))
		one.pt <- function(rg.pt) {
			wgts <- kernel( (x-rg.pt)/bw )/bw
			wgts %*% y /sum(wgts)
		}
		yp <- sapply(x.points, one.pt)
	}
	list(x=x.points, y=yp)
}