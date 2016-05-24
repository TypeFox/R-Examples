# computing probability of all links of the graph
plinks = function( bdgraph.obj, round = 3 )
{
	p_links = bdgraph.obj $ p_links

	if( !is.null( bdgraph.obj $ sample_graphs ) )
	{
		sample_graphs <- bdgraph.obj $ sample_graphs
		graph_weights <- bdgraph.obj $ graph_weights
		p             <- nrow( bdgraph.obj $ last_graph )
		pvec          <- c( rep( 0, p * ( p - 1 ) / 2 ) )
	   
		for ( i in 1 : length( sample_graphs ) )
		{
			inp       <- which( unlist( strsplit( as.character( sample_graphs[i] ), "" ) ) == 1 )
			pvec[inp] <- pvec[inp] + graph_weights[i]
		}
		
		dimlab  <- colnames( bdgraph.obj $ last_graph ) # lastG
		p_links <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )
		p_links[ upper.tri(p_links) ] <- pvec / sum( graph_weights )
	}
	
	return( Matrix( round( p_links, round ) ) )
}
      
