plugin.select <- function(x, y, drv = 0L, degree = drv+1L) 
{	
	n <- nrow(y)
	N <- ncol(y)
	if (length(x) != N) 
		stop("length(x) and nrow(y) must match")
	if (!(drv %in% c(0,1))) 
		stop("drv must be 0 or 1")
	if (!((degree-drv) %in% c(0,1))) 
		stop("(degree - drv) must be 0 or 1")
	spacing <- diff(x)
	if (any(spacing < 0) || !isTRUE(all.equal(min(spacing),max(spacing))))
		stop("x must be a uniform grid") 
	range.x <- x[N] - x[1]
	int.K1 	<- 1.1283
	int.K2 	<- 0.5641
	y.mean  <- colMeans(y, na.rm=TRUE)
	I.alpha <- mean(apply(y[,-1]-y[,-N], 1, crossprod))
	h.cv 	<- cv.select(x, y, degree)
	mu.cv 	<- locpoly(x, y.mean, degree = degree, bandwidth = h.cv, gridsize = N)$y
	Dmu.cv  <- diff(mu.cv, differences = drv + 2) * (N-1)^(drv+2)
	L2 		<- sum(Dmu.cv^2) / (N-1)
	if (drv == 0) h.plug <- max((I.alpha * int.K1 / 2 / L2 / n)^(1/3), .5/(N-1))
	if (drv == 1) h.plug <- max((I.alpha * int.K2 / 2 / L2 / n)^(1/5), .5/(N-1))	
	return(h.plug * range.x)
}
