
clv.Scatt <- function( data, clust, dist="euclidean")
{
	data = data.validity(data, "data")
	clust = cls.id.vect.validity(clust, "clust")
	dist = dist.validity(dist)
	
	if(dim(data)[1] != length(clust))
		stop("Bad input data: number of 'data' objects do not agree with length of vector 'clust'.")
	
	clust_num = as.integer(max(clust))
	scatt <- .Call( "Scatt",
					data,
					clust,
					clust_num,
					dist,
					PACKAGE="clv"
				  )

	class(scatt) = "scatt.obj"

	return(scatt)
}

clv.Dis <- function(cluster.center)
{
	cluster.center = data.validity(cluster.center, "cluster.center")

	dis <- .Call( "Dis",
				  cluster.center,
				  PACKAGE="clv"
				)

	return(dis)
}

clv.DensBw <- function( data, clust, scatt.obj, dist="euclidean" )
{
	data = data.validity(data, "data")
	clust = cls.id.vect.validity(clust, "clust")
	dist = dist.validity(dist)
	
	if( class(scatt.obj) != "scatt.obj" )
		stop("Bad input data: 'scatt.obj' is not a result of 'clv.Scatt' function.")
	
	if(dim(data)[1] != length(clust))
		stop("Bad input data: number of 'data' objects do not agree with length of vector 'clust'.")

	dens <- .Call( "Dens_bw",
					data,
					clust,
					scatt.obj$cluster.center,
					scatt.obj$stdev,
					dist,
					PACKAGE="clv"
				  )

	return(dens)
}

clv.SD <- function( scatt, dis, alfa )
{
	return( alfa*scatt + dis )	
}

clv.SDbw <- function( scatt, dens )
{
	return( scatt + dens )	
}
