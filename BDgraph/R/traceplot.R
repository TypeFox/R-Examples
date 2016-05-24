# plot of graph size to check the convergency of BDMCMC algorithm
traceplot = function( bdgraph.obj, acf = FALSE, pacf = FALSE, main = NULL, ... )
{
	if( !is.null( bdgraph.obj $ p_links ) ) stop( "Function needs output of 'bdgraph' with option save.all = TRUE" )  

	sample_graphs     = bdgraph.obj $ sample_graphs
    all_graphs        = bdgraph.obj $ all_graphs
	graph_weights     = bdgraph.obj $ graph_weights
	sizesample_graphs = sapply( sample_graphs, function(x) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )  
	sizeall_graphs    = sizesample_graphs[ all_graphs ]
	which_G_max      = which( max( graph_weights ) == graph_weights )
	size_selected_g  = sizeall_graphs[ which_G_max ] 

	if( is.null( main ) ) main = "Trace of graph size"
	
	x_vec = 1 : length( all_graphs )
	
	if ( acf == FALSE & pacf == FALSE )
	{
		plot( x = x_vec, sizeall_graphs, type = "l", main = main, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_g, col = "red" )	   
	}
	
	if ( acf == TRUE & pacf == TRUE )
	{
		op = par( mfrow = c( 2, 2 ), pty = "s" )  
		plot( x = x_vec, sizeall_graphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_g, col = "red" )	  
		acf( sizeall_graphs,  main = "ACF for graph size" )
		pacf( sizeall_graphs, main = "PACF for graph size" )
		par( op )
	}
	
	if ( acf == TRUE & pacf == FALSE )
	{
		op <- par( mfrow = c( 1, 2 ), pty = "s" ) 
		plot( x = x_vec, sizeall_graphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_g, col = "red" )	  
		acf( sizeall_graphs, main = "ACF for graph size" )
		par( op )
	}
	
	if ( acf == FALSE & pacf == TRUE )
	{
		op <- par( mfrow = c( 1, 2 ), pty = "s" ) 
		plot( x = x_vec, sizeall_graphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		abline( h = size_selected_g, col = "red" )	  
		pacf( sizeall_graphs, main = "PAIC for graph size" )
		par(op)
	}		
}  
      
