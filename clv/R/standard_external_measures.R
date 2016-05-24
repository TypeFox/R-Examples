

# faster alghorithm

std.ext <- function(clust1, clust2)
{
	clust1 = cls.id.vect.validity(clust1, "clust1")
	clust2 = cls.id.vect.validity(clust2, "clust2")

	if( length(clust1) != length(clust2) )
		stop("Bad input data: both vectors should have the same length.")

	cls.num1 = as.integer(max(clust1))
	cls.num2 = as.integer(max(clust2))
	
	cnf.mx = .Call("confusionMatrix",
					clust1,
					clust2,
					c(cls.num1,cls.num2),
					PACKAGE="clv"
				   )

	result = .Call("standardExternalMeasures",
					cnf.mx,
					PACKAGE="clv")
	return( result )
}

# slower alghorithm

ext.measures.slow <- function(clust1, clust2)
{
	clust1 = cls.id.vect.validity(clust1, "clust1")
	clust2 = cls.id.vect.validity(clust2, "clust2" )

	if( length(clust1) != length(clust2) )
		stop("Bad input data: both vectors should have the same length.")

	result = .Call("standardExternalMeasuresSlow",
					clust1,
					clust2,
					PACKAGE="clv"
				   )
	return( result )
}

# function convert check if it is a vector or list produced by 'ext.measures' function
# and convert it to a vector

clv_conv <- function(external.ind)
{
	result = NA
	if( is.vector(external.ind) & is.numeric(external.ind)  & length(external.ind) == 4 )
	{
		result = external.ind
	}
	else
	{
		if( is.list(external.ind) )
		{
			result = as.integer( c(external.ind$SS, external.ind$SD, external.ind$DS, external.ind$DD) )
			if( length(result) != 4 )
				stop( "Bad input data: 'external.ind' is not list produced by 'ext.measures' function.")
		}
		else stop( "Bad input data: 'external.ind' is neither integer vector nor list type." )
	}
	return(result)
}

# input data in functions presented below 
# is an object which is produced by "ext.measures"
# or a integer vector where all four and only four values are > 0 

clv.Rand <- function(external.ind)
{
	v = clv_conv(external.ind)
	result = (v[1]+v[4])/sum(v)
	return(result)
}

clv.Jaccard <- function(external.ind)
{
	v = clv_conv(external.ind)
	result = v[1]/(v[1] + v[2] + v[3])
	return(result)
}

clv.Folkes.Mallows <- function(external.ind)
{
	v = clv_conv(external.ind)
	m1 = v[1]/( v[1] + v[2] )
	m2 = v[1]/( v[1] + v[3] )
	result = sqrt( m1*m2 )
	return(result)
}

clv.Russel.Rao <- function(external.ind)
{
	v = clv_conv(external.ind)
	result = v[1]/sum(v)
	return(result)
}

clv.Phi <- function(external.ind)
{
	v = clv_conv(external.ind)
	result = (v[1]*v[4] - v[2]*v[3])/sqrt((v[1]+v[2])*(v[1]+v[3])*(v[2]+v[4])*(v[3]+v[4])) 
	return(result)
}

