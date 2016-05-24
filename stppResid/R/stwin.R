stwin <- function(xcoord = c(0, 1), ycoord = c(0, 1), tcoord = c(0, 1))
{
	if (!is.vector(xcoord) || length(xcoord) != 2 || xcoord[2] < xcoord[1])
		stop("xcoord must be a vector of minimum and maximum x-coordinates (min, max)")
	if (!is.vector(ycoord) || length(ycoord) != 2 || ycoord[2] < ycoord[1])
		stop("ycoord must be a vector of minimum and maximum y-coordinates (min, max)")
	if (!is.vector(tcoord) || length(tcoord) != 2 || tcoord[2] < tcoord[1])
		stop("tcoord must be a vector of minimum and maximum t-coordinates (min, max)")
	stw <- list(xcoord = xcoord, ycoord = ycoord, tcoord = tcoord)
	class(stw) <- "stwin"
	return(stw)	
}