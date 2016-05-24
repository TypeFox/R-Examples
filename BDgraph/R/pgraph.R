# computing the probability of all the possible graphs or one specific graph 
pgraph = function( bdgraph.obj, number.g = 4, adj_g = NULL )
{
	if( !is.null( bdgraph.obj $ p_links ) ) stop( "It needs object of 'bdgraph' with option save.all = TRUE" ) 
	
	sample_graphs  = bdgraph.obj $ sample_graphs
	graph_weights  = bdgraph.obj $ graph_weights
	sort_gWeights = sort( graph_weights, decreasing = TRUE )
  
	if( is.null( adj_g ) )
	{
		p      <- nrow( bdgraph.obj $ last_graph )
		list_g <- list()
		vec_g  <- c( rep( 0, p * ( p - 1 ) / 2 ) )  

		for ( i in 1 : number.g )
		{
			vec_g <- 0 * vec_g
			indG_i = sample_graphs[ which( graph_weights == sort_gWeights[i] ) ]
			vec_g[ which( unlist( strsplit( as.character( indG_i ), "" ) ) == 1 ) ] <- 1
			list_g[[i]] <- matrix( 0, p, p )
			list_g[[i]][ upper.tri( list_g[[i]] ) ] <- vec_g
			list_g[[i]] <- Matrix( list_g[[i]], sparse = TRUE )
		}

		return( list( selected_g = list_g, prob_g = sort_gWeights[1 : number.g] / sum( graph_weights ) ) )
		
	} 
	else 
	{
		if ( class(adj_g) == "sim" ) G <- as.matrix( adj_g $ G )

		indG   = paste( G[upper.tri(G)], collapse = '' )
		wh     = which( sample_graphs == indG )
		prob_g = ifelse( length(wh) == 0, 0, graph_weights[wh] / sum( graph_weights ) )

		return( prob_g )
	}
}
    
