

wcls.matrix <- function(data, clust, cluster.center)
{
	data = data.validity(data, "data")
	clust = cls.id.vect.validity(clust, "clust")
	cluster.center = data.validity(cluster.center, "cluster.center")

	if( dim(data)[1] != length(clust) )
		stop("Bad input data: number of 'data' objects is not the same as lenght of 'cluster' vector.")
	if( dim(data)[2] != dim(cluster.center)[2] )
		stop("Bad input data: dimension of 'data' objects is not the same as dimension of 'cluster' centers.")

	result = .Call("whithinClusterScatterMatrix", 
					data, 
					clust,
					cluster.center,
					PACKAGE="clv")
	return( result )
}

bcls.matrix <- function( cluster.center, cluster.size, mean )
{
	cluster.center = data.validity(cluster.center, "cluster.center")
	if( !is.vector(mean) || !is.numeric(mean) )
		stop("Bad input data: 'mean' is not a numeric vector type.")
	if( !is.vector(cluster.size) )
		stop("Bad input data: 'cluster.size' is not a vector type.")
	if( dim(cluster.center)[1] != length(cluster.size) )
		stop("Bad input data: number of 'cluster.center' objects is not the same as lenght of 'cluster.size' vector.")
	if( dim(cluster.center)[2] != length(mean) )
		stop("Bad input data: dimension of 'data' objects is not the same as dimension of 'cluster' centers.")

	result = .Call("betweenClusterScatterMatrix", 
					cluster.center, 
					as.integer(cluster.size),
					mean,
					PACKAGE="clv")
	return( result )
}
