add.stpoints <- function(X, points)
{
	if(!is.stpp(X))
		stop("X must be an object of type stpp")
	if(!is.vector(points) && !is.data.frame(points) && !is.matrix(points) && !is.list(points))
		stop("points must be a vector, data frame, matrix, or list")
	if(is.vector(points) || is.list(points)) {
		if(length(points) != 3)
			stop("points vector or list must be of length 3 (x, y, t)")
		if(is.list(points)){
			if((length(points[[1]]) != length(points[[2]])) || (length(points[[1]]) != length(points[[3]])) || (length(points[[2]]) != length(points[[3]])))
				stop("elements of points list must be of same length")
			Y <- stpp(c(X$x, points[[1]]), c(X$y, points[[2]]), c(X$t, points[[3]]), stwin(X$xcoord, X$ycoord, X$tcoord))
		} else
			Y <- stpp(c(X$x, points[1]), c(X$y, points[2]), c(X$t, points[3]), stwin(X$xcoord, X$ycoord, X$tcoord))
	}
	if(is.data.frame(points) || is.matrix(points)) {
		if(ncol(points) < 3)
			stop("points data frame must contain 3 columns [x, y, t]")
		Y <- stpp(c(X$x, points[,1]), c(X$y, points[,2]), c(X$t, points[,3]), stwin(X$xcoord, X$ycoord, X$tcoord))
	}
	return(Y)
}