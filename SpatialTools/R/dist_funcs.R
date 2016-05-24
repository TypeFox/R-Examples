#Fast C function for Euclidean distance of a matrix of coordinates
dist1 <- function(coords)
{
	#check validity fo coords argument
	if(!is.matrix(coords) || !is.numeric(coords))
	{
		stop("coords must be a numeric matrix")
	}
	
	nr <- nrow(coords)
	matrix(.C("dist1", x = as.double(coords), nc = as.integer(ncol(coords)), 
		nr = as.integer(nr), 
		d = as.double(numeric(nr^2)), PACKAGE = "SpatialTools")$d, nrow = nr)
}

#Fast C function for Euclidian distance between two pairs of matrices.
#If coords1 is of dimension n1 x r, and coords2 of dimension n2 x r,
#then the resulting distance matrix is of size n1 x n2
dist2 <- function(coords1, coords2)
{
	#check validity fo coords1 and coords2 arguments
	
	if(!is.matrix(coords1) || !is.numeric(coords1))
	{
		stop("coords1 must be a numeric matrix")
	}
	if(!is.matrix(coords2) || !is.numeric(coords2))
	{
		stop("coords2 must be a numeric matrix")
	}
	if(ncol(coords1) != ncol(coords2))
	{
		stop("coords1 and coords2 must have the same number of columns")
	}

	nr1 <- nrow(coords1)
	nr2 <- nrow(coords2)
	matrix(.C("dist2", x1 = as.double(coords1), x2 = as.double(coords2), 
		nc = as.integer(ncol(coords1)), nr1 = as.integer(nr1), 
		nr2 = as.integer(nr2), 
		d = as.double(numeric(nr1*nr2)), PACKAGE = "SpatialTools")$d, nrow = nr1)
}

#Function to determine coincident locations
coincident <- function(coords1, coords2)
{
	#check arguments of coincident locations
	if(!is.matrix(coords1) || !is.numeric(coords1))
	{
		stop("coords1 must be a numeric matrix")
	}
	if(!is.matrix(coords2) || !is.numeric(coords2))
	{
		stop("coords2 must be a numeric matrix")
	}
	if(ncol(coords1) != ncol(coords2))
	{
		stop("coords1 and coords2 must have the same number of columns")
	}
	
	.Call( "coincident_cpp", coords1 = coords1, coords2 = coords2, 
		eps = .Machine$double.eps, PACKAGE = "SpatialTools")
}