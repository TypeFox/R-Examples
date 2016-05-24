
cls.set.section <- function (clust1,clust2)
{
	clust1 = data.validity.int(clust1, "clust1")
	clust2 = data.validity.int(clust2, "clust2")
 
	if( dim(clust1)[2] != 2 || dim(clust2)[2] != 2 )
		stop("Bad input data: each matrix should have two columns.")

	result <- .Call( "clv_clusteredSetsSection" ,
					as.integer(clust1),
					as.integer(clust2),
					as.integer(c(dim(clust1)[1],dim(clust2)[1])),
					PACKAGE="clv" )

	return(result)
}

# distance between two partitioning of the same data set

dot.product <- function (clust1,clust2) 
{
	clust1 = cls.id.vect.validity(clust1,"clust1")
	clust2 = cls.id.vect.validity(clust2,"clust2")
 
	if( length(clust1) != length(clust2) )
		stop("Bad input data: both vectors should have the same length.")

	result <- .Call( "clv_dotProduct" ,
					clust1,
					clust2,
					PACKAGE="clv" )

	if( is.infinite(result) ) result = NaN

	return(result)
}
