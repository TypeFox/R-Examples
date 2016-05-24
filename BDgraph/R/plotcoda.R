# To check the convergency of the BDMCMC algorithm
plotcoda = function( bdgraph.obj, thin = NULL, control = TRUE, main = NULL, ... )
{
	if( !is.null( bdgraph.obj $ p_links ) ) stop( "It needs object of 'bdgraph' with option save.all = TRUE" ) 
	
	if( is.null( thin ) ) thin = ceiling( length( bdgraph.obj $ all_graphs ) / 1000 )

	sample_graphs   = bdgraph.obj $ sample_graphs
	p               = nrow( bdgraph.obj $ last_graph )
	qp              = p * ( p - 1 ) / 2 
	all_weights     = bdgraph.obj $ all_weights
	all_graphs      = bdgraph.obj $ all_graphs

	allG_new        = all_graphs[ c( thin * ( 1 : floor( length( all_graphs ) / thin ) ) ) ]
	all_weights_new = all_weights[ c( thin * ( 1 : floor( length( all_weights ) / thin ) ) ) ]
	length_allG_new = length( allG_new )
	result          = matrix( 0, qp, length_allG_new )
	vec_result      = 0 * result[ , 1]

	for ( g in 1 : length_allG_new )
	{
		mes = paste( c( "Calculation ... in progress : ", floor( 100 * g / length_allG_new ), "%" ), collapse = "" )
		cat(mes, "\r")
		flush.console()	

		which_edge             = which( unlist( strsplit( as.character( sample_graphs[ allG_new[g] ] ), "" ) ) == 1 )
		vec_result[which_edge] = vec_result[which_edge] + all_weights_new[g]
		result[ ,g]            = vec_result / sum( all_weights_new[ c( 1 : g ) ] )    	 
	}

	if ( control )
		if ( p > 15 )
		{
			randomLinks = sample( x = 1:qp, size = ( qp - 100 ), replace = FALSE )
			result[ randomLinks, ] = 0
		}
	
	mes = paste( c( "Calculation ... done.                        " ), collapse = "" )
	cat(mes, "\r")
	cat("\n")
	flush.console()

	matplot( x = thin * ( 1 : length_allG_new ), y = t( result ), type = "l", lty = 1, col = 1,
		  xlab = "Iteration", ylab = "Posterior link probability", cex.lab = 1.3, cex.axis = 1.2 )
		  
	if ( is.null( main ) ) main = "Trace of the Posterior Probabilities of the Links"
	title( main = main, cex.main = 1.5 )
}
    
