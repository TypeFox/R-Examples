# To compare the result
ROC = function ( G, est ) 
{
	tp   = sum( ( G != 0 ) * ( est != 0 ) )
	fp   = sum( ( G == 0 ) * ( est != 0 ) )
	fn   = sum( ( G != 0 ) * ( est == 0 ) )
	tn_M = ( G == 0 ) * ( est == 0 )
	tn   = sum( tn_M[ upper.tri( tn_M == 1 ) ] )
	
	# Precision is the probability that a randomly selected link is relevant 
	Precision <- tp / ( tp + fp ) 
	
	# Recall is the probability that a randomly selected relevant link 	
	Recall   <- tp / ( tp + fn )               # also called TPR
	FPR      <- fp / ( fp + tn )               # False positive rate
	Accuracy <- ( tp + tn ) / ( tp + tn + fp + fn )
	
	# harmonic mean of precision and recall, called F-measure or balanced F-score
	F1score <- ( 2 * tp ) / ( 2 * tp + fp + fn )
	
	return( c( tp, tn, fp, fn, Recall, FPR, Accuracy, F1score, Precision ) )
}
   
# To compare the result according to the true graph
compare = function ( sim.obj, bdgraph.obj, bdgraph.obj2 = NULL, bdgraph.obj3 = NULL, colnames = NULL, vis = FALSE ) 
{
	if( is.matrix( sim.obj ) )      G    = sim.obj
	if( is.matrix( bdgraph.obj ) )  est  = bdgraph.obj
	if( is.matrix( bdgraph.obj2 ) ) est2 = bdgraph.obj2
	if( is.matrix( bdgraph.obj3 ) ) est3 = bdgraph.obj3
	
	if( class(sim.obj)      == "sim" )      G    <- sim.obj $ G 
#~ 	if( class(sim.obj)      == "bdgraph" ){ es   <- select( sim.obj ); G <- est $ G ; est <- es }
	if( class(bdgraph.obj)  == "bdgraph" )  est  <- select( bdgraph.obj ) 
	if( class(bdgraph.obj2) == "bdgraph" )  est2 <- select( bdgraph.obj2 ) 
	if( class(bdgraph.obj3) == "bdgraph" )  est3 <- select( bdgraph.obj3 ) 
	if( class(bdgraph.obj)  == "select"  )  est  <- bdgraph.obj $ refit
	if( class(bdgraph.obj2) == "select"  )  est2 <- bdgraph.obj2 $ refit
	if( class(bdgraph.obj3) == "select"  )  est3 <- bdgraph.obj3 $ refit
   
	G   = as.matrix(G)        # G is the adjacency matrix of true graph 
	est = as.matrix(est)      # est is the adjacency matrix of estimated graph 
	p   = nrow( G )
	G[ lower.tri( G, diag = TRUE ) ]     = 0
	est[ lower.tri( est, diag = TRUE ) ] = 0
   	   
	if( is.null( bdgraph.obj2 ) & is.null( bdgraph.obj3 ) )
	{
		result = matrix( 1, 9, 2 )
		result[ , 2 ] = ROC( G = G, est = est )
		if( is.null( colnames ) ) colnames <- c( "True graph", "estimate" )
	}

	if( !is.null( bdgraph.obj2 ) & is.null( bdgraph.obj3 ) )
	{
		est2 = as.matrix( est2 )       
		est2[ lower.tri( est2, diag = TRUE ) ] = 0

		result = matrix( 1, 9, 3 )
		result[ , 2 ] = ROC( G = G, est = est )
		result[ , 3 ] = ROC( G = G, est = est2 )
		if( is.null( colnames ) ) colnames <- c( "True graph", "estimate", "estimate2" )
	} 
	
	if( !is.null( bdgraph.obj2 ) & !is.null( bdgraph.obj3 ) )
	{
		est2 = as.matrix( est2 )       
		est2[ lower.tri( est2, diag = TRUE ) ] = 0
		est3 = as.matrix( est3 )       
		est3[ lower.tri( est3, diag = TRUE ) ] = 0

		result = matrix( 1, 9, 4 )
		result[ , 2 ] = ROC( G = G, est = est )
		result[ , 3 ] = ROC( G = G, est = est2 )
		result[ , 4 ] = ROC( G = G, est = est3 )
		if( is.null( colnames ) ) colnames <- c( "True graph", "estimate", "estimate2", "estimate3" )
	} 
	
   if( vis )
   {
		G_igraph   <- graph.adjacency( G,   mode = "undirected", diag = FALSE )
		est_igraph <- graph.adjacency( est, mode = "undirected", diag = FALSE )
		if ( p < 20 ) sizev = 15 else sizev = 2

		if( is.null( bdgraph.obj2 ) )
		{
			op <- par( mfrow = c( 1, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )
			plot.igraph( G_igraph,   layout = layout.circle, main = colnames[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est_igraph, layout = layout.circle, main = colnames[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		}
		 
		if( !is.null( bdgraph.obj2 ) & is.null( bdgraph.obj3 ) )
		{
			op   <- par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )
			est2_igraph <- graph.adjacency( as.matrix(est2), mode = "undirected", diag = FALSE )
			plot.igraph( G_igraph,    layout = layout.circle, main = colnames[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est_igraph,  layout = layout.circle, main = colnames[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est2_igraph, layout = layout.circle, main = colnames[3], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		if( !is.null( bdgraph.obj2 ) & !is.null( bdgraph.obj3 ) )
		{
			op   <- par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )
			est2_igraph <- graph.adjacency( as.matrix(est2), mode = "undirected", diag = FALSE ) 
			est3_igraph <- graph.adjacency( as.matrix(est3), mode = "undirected", diag = FALSE ) 
			plot.igraph( G_igraph,    layout = layout.circle, main = colnames[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est_igraph,  layout = layout.circle, main = colnames[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
			plot.igraph( est2_igraph, layout = layout.circle, main = colnames[3], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
			plot.igraph( est3_igraph, layout = layout.circle, main = colnames[4], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		par(op)
   }
   
	result[ c( 3, 4, 6 ),1 ] = 0
	result[ 1, 1 ] = sum( G )
	result[ 2, 1 ] = p * ( p - 1 ) / 2 - result[ 1, 1 ]  

	result[ is.na( result ) ] = 0

	colnames( result ) <- colnames
	rownames( result ) <- c("true positive", "true negative", "false positive", "false negative", 
                "true positive rate", "false positive rate", "accuracy", "balanced F-score", "positive predictive")
				
	return( round( result, 3 ) )
}
         
