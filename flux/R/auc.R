"auc" <- function(x, y, thresh = NULL, dens = 100, sort.x = TRUE){
	##calculates the area under curve (integral) for any given set of x- and y-coordinates
	# check for and omit na entries in x or y
	x <- x[!is.na(x)]
	y <- y[!is.na(x)]
	x <- x[!is.na(y)]
	y <- y[!is.na(y)]
	##order after x (default)
	if(sort.x){
		ord <- order(x)
		x <- x[ord]
		y <- y[ord]
	}
	# increase density of points
	idx = 2:length(x)
	x <- as.vector(apply(cbind(x[idx-1], x[idx]), 1, function(x) seq(x[1], x[2], length.out=dens)))
	y <- as.vector(apply(cbind(y[idx-1], y[idx]), 1, function(x) seq(x[1], x[2], length.out=dens)))
	# acknowledge threshold
	if(!is.null(thresh)){
		y.0 <- y <= thresh
		y[y.0] <- thresh
	}
	# calculate area under curve
	idx = 2:length(x)
	integral <- as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
    integral
}