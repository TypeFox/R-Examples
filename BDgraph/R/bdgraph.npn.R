# non-parametric transfer function for non-normal data
bdgraph.npn = function( data, npn = "shrinkage", npn.thresh = NULL )
{
    if( !is.matrix( data ) & !is.data.frame( data ) ) stop( "Data should be a matrix or dataframe" )	
    if( is.data.frame( data ) ) data <- data.matrix( data )
    if( any( is.na( data ) ) ) stop( "Data should contain no missing data" ) 
	
	n <- nrow( data )
	
  	# shrinkage transfer
	if( npn == "shrinkage" )
	{
		data <- qnorm( apply( data, 2, rank ) / ( n + 1 ) )
		data <- data / sd( data[ , 1] )
	}
	
	# truncation transfer
	if( npn == "truncation" )
	{
		if( is.null( npn.thresh ) ) npn.thresh <- 0.25 * ( n ^ - 0.25 ) * ( pi * log(n) ) ^ - 0.5
		
		data <- qnorm( pmin( pmax( apply( data, 2, rank ) / n, npn.thresh ), 1 - npn.thresh ) )
    	data <- data / sd( data[ , 1] )
	}

	if( npn == "skeptic" ) data <- 2 * sin( pi / 6 * cor( data, method = "spearman" ) )
	
	return( data )
}
