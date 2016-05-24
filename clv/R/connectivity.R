
connectivity <- function ( data, clust, neighbour.num, dist="euclidean" )
{
	data = data.validity(data, "data")
	clust = cls.id.vect.validity(clust, "clust")
	dist = dist.validity(dist)	
 
	if( is.vector(neighbour.num) == FALSE ) 
		stop("Bad usage: input 'neighbour.num' should be vector type.")
	
	if(dim(data)[1] != length(clust))
		stop("Bad input data: number of 'data' objects do not agree with length of vector 'clust'.")

	if( !is.numeric(neighbour.num) )
		stop("Bad input data: 'neighbour.num' is not numeric type.")
	if( neighbour.num[1] < 1 | neighbour.num[1] > length(clust) )
		stop("Bad input data: parameter 'neighbour.num' should be greater than 0 and less than number of analized objects.")

	result <- .Call("connectivity",
					data,
					clust,
					as.integer(neighbour.num),
					dist,
					PACKAGE="clv"
					)
	return(result)
}

connectivity.diss.mx <- function ( diss.mx, clust, neighbour.num )
{
	diss.mx = data.validity(diss.mx, "diss.mx")
	clust = cls.id.vect.validity(clust, "clust")
 
	if( is.vector(neighbour.num) == FALSE ) 
		stop("Bad usage: input 'neighbour.num' should be vector type.")
	
	if(dim(diss.mx)[1] != dim(diss.mx)[2])
		stop("Bad input data: 'diss.mx' should be a square (symetric) matrix.")

	if(dim(diss.mx)[1] != length(clust))
		stop("Bad input data: number of 'diss.mx' objects do not agree with length of vector 'clust'.")

	if( !is.numeric(neighbour.num) )
		stop("Bad input data: 'neighbour.num' is not numeric type.")
	if( neighbour.num[1] < 1 | neighbour.num[1] > length(clust) )
		stop("Bad input data: parameter 'neighbour.num' should be greater than 0 and less than number of analized objects.")

	result <- .Call("connectivityDissMx",
					diss.mx,
					clust,
					as.integer(neighbour.num),
					PACKAGE="clv"
					)
	return(result)
}
