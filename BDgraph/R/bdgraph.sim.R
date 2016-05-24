# Data generator according to the graph structure
bdgraph.sim = function( n = 2, p = 10, type = "Gaussian", graph = "random", 
						prob = 0.2, size = NULL, mean = 0, class = NULL, 
						cut = 4, b = 3, D = diag(p), K = NULL, sigma = NULL, 
						vis = FALSE )
{
    if ( is.matrix(K) )  graph <- "fixed"
    
    if ( type == "normal" )     type = "Gaussian"
    if ( type == "non-normal" ) type = "non-Gaussian"
    
    if ( is.matrix(graph) )
	{
		G     <- graph
	    graph <- "fixed"
    } 
	
	# build the graph structure
	if ( graph == "random" )
	{
		G <- matrix( 0, p, p )

		if ( is.null(size) )
		{
			if ( prob < 0 | prob > 1 ) stop( "'prob' should be between zero and one" )
			
			G[upper.tri(G)] <- rbinom( p * ( p - 1 ) / 2, 1, prob )
		} 
		else 
		{
			if ( size < 0 | size > p * ( p - 1 ) / 2 )  stop( "Graph size should be between zero and p*(p-1)/2" )
			
			smp <- sample( 1 : ( p * ( p - 1 ) / 2 ), size, replace = FALSE )
			G[upper.tri(G)][smp] <- 1
		}
		
		G <- G + t(G)
	}
	
	if ( graph == "cluster" )
	{
		# partition variables
		if ( is.null(class) )
		{ 
			if ( !is.null(size) )
			{
				class <- length(size)
			} 
			else 
			{
				class <- max( 2, ceiling( p / 20 ) )
			}
		}

		g.large <- p %% class
		g.small <- class - g.large
		n.small <- floor( p / class )
		n.large <- n.small + 1
		vp      <- c( rep( n.small, g.small ), rep( n.large, g.large ) )
		 
		G       <- matrix( 0, p, p )
		 
		if ( is.null(size) )
		{
			if ( prob < 0 | prob > 1 ) stop( "'prob' should be between zero and one" )

			for ( i in 1 : class )
			{
				tmp <- if ( i == 1 ) (1 : vp[1]) else ( ( sum( vp[1 : (i-1)] ) + 1 ) : sum( vp[1:i] ) )
				gg                <- matrix( 0, vp[i], vp[i] )
				gg[upper.tri(gg)] <- rbinom( vp[i] * (vp[i] - 1) / 2, 1, prob )
				G[tmp, tmp]       <- gg
			}
		} 
		else 
		{
			if ( class != length(size) )  stop( "Number of graph sizes is not match with number of clusters" )
			if ( sum(size) < 0 | sum(size) > p * (p - 1) / 2 )   stop( "Total graph sizes should be between zero and p*(p-1)/2" )

			for (i in 1 : class)
			{
				tmp <- if ( i == 1 ) ( 1 : vp[1] ) else ( ( sum( vp[1 : (i-1)] ) + 1 ) : sum( vp[1:i] ) )
				gg  <- matrix( 0, vp[i], vp[i] )
				smp <- sample( 1 : ( vp[i] * (vp[i] - 1) / 2 ), size[i], replace = FALSE )
				gg[upper.tri(gg)][smp] <- 1
				G[tmp, tmp]            <- gg
			}
		}
		
		G <- G + t(G)	   
	}

	if( graph == "hub" )
	{
		# partition variables
		if ( is.null( class ) ) class <- ceiling( p / 20 ) 

		# partition variables into groups
		g.large <- p %% class
		g.small <- class - g.large
		n.small <- floor( p / class )
		n.large <- n.small + 1
		g.list  <- c( rep( n.small, g.small ), rep( n.large, g.large ) )
		g.ind   <- rep( c( 1:class ), g.list )
		
		G <- matrix(0, p, p)
		
		for ( i in 1:class )
		{
			tmp           <- which( g.ind == i )
			G[tmp[1],tmp] <- 1
			G[tmp,tmp[1]] <- 1
		}
	}
	
	if ( graph == "circle" )
	{
	    G       <- toeplitz( c( 0, 1, rep( 0, p - 2 ) ) )
        G[1, p] <- 1
		G[p, 1] <- 1
	}

	if ( graph == "scale-free" )
	{
#		data_sim <- huge.generator( n = 2, d = 5, graph = "scale-free" )
		G = matrix( 0, p, p )
		# scaleFree( int *G, int *p )
		resultGraph = .C( "scaleFree", G = as.integer(G), as.integer(p), PACKAGE = "BDgraph" )
		G           = matrix( resultGraph $ G, p, p ) 
	}
	
    if ( !is.null(sigma) ) K <- solve(sigma)   

    if ( is.matrix(K) )
    { 
		G     <- 1 * ( abs(K) > 0.02 )
		if( is.null( sigma ) ) sigma <- solve(K)	
    } 
    else 
    {
		Ti        = chol( solve(D) )
		diag(G)   = 0
		K         = matrix( 0, p, p )
		threshold = 1e-8
		
		result = .C( "rgwish", as.integer(G), as.double(Ti), K = as.double(K), 
		             as.integer(b), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
		
		K     = matrix ( result $ K, p, p ) 		
		sigma = solve( K )
	}
	
	diag(G) <- 0
	p       <- nrow( G )
	
	# generate multivariate normal data
	if ( typeof(mean) == "double" ) mean <- rep( mean, p )
	R <- chol( sigma )
	z <- matrix( rnorm( p * n ), p, n )
	d <- t(R) %*% z + mean
	d <- t(d)

	if ( type == "mixed" )
	{
		# generating mixed data which are 'count', 'ordinal', 'non-Gaussian', 
		# 'binary', and 'Gaussian', respectively.
		ps = floor( p / 5 )
		
		# generating count data
		col_number     <- c( 1:ps )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qpois( p = prob, lambda = 10 )

		# generating ordinal data
		col_number     <- c( ( ps + 1 ):( 2 * ps ) )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qpois( p = prob, lambda = 2 )

		# generating non-Guassian data
		col_number     <- c( (2 * ps + 1):(3 * ps) )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qexp( p = prob, rate = 10 )

		# for binary data
		col_number     <- c( (3 * ps + 1):(4 * ps) )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qbinom(p = prob, size = 1, prob = 0.5 )
	}

	if ( type == "non-Gaussian" )
	{
		# generating multivariate continuous non-Gaussian data  
		prob <- pnorm( d )
		d    <- qexp( p = prob, rate = 10 )
	}

	if ( type == "discrete" )
	{
		runif_m   <- matrix( runif( cut * p ), nrow = p, ncol = cut )   
		marginals <- apply( runif_m, 1, function(x) { qnorm( cumsum( x / sum(x) )[-length(x)] ) } )
			 
		for ( j in 1:p )
		{
			breaks <- c( min( d[,j] ) - 1, marginals[,j], max(d[,j]) + 1 )  
			d[,j]  <- as.integer( cut( d[,j], breaks = breaks, right = FALSE ) )
		}	
	}
	
	# graph visualization
	if ( vis )
	{
		graphG <- graph.adjacency( G, mode = "undirected", diag = FALSE )
		
		if ( p < 20 ) size = 10 else size = 2
		plot.igraph( graphG, layout = layout.circle, main = "Graph structure", 
		             vertex.color = "white", vertex.size = size, vertex.label.color = 'black' )
	}
	
	dimlab     <- as.character( 1 : p )
	simulation <- list( G = Matrix( G, sparse = TRUE, dimnames = list( dimlab, dimlab ) ), 
	                    data = d, sigma = sigma, K = K, graph = graph, type = type )
	
	class( simulation ) <- "sim"
	return( simulation )
}
    
# Print function for simulation data
print.sim = function( x, ... )
{
	p <- ncol( x $ sigma )

	cat( paste( "  Data generated by bdgraph.sim"                           ), fill = TRUE )
	cat( paste( "  Data type       =", x $ type                             ), fill = TRUE )
	cat( paste( "  Sample size     =", nrow(x $ data)                       ), fill = TRUE )
	cat( paste( "  Graph type      =", x $ graph                            ), fill = TRUE )
	cat( paste( "  Number of nodes =", p                                    ), fill = TRUE )
	cat( paste( "  Graph size      =", sum(x $ G) / 2                       ), fill = TRUE )
	cat( paste( "  Sparsity        =", round(sum(x $ G) / (p * (p - 1)), 4) ), fill = TRUE )
}
    
# plot for class "sim" from bdgraph.sim function
plot.sim = function( x, main = NULL, layout = layout.circle, ... )
{
    true_graph = as.matrix( x $ G )
    if ( is.null( main ) ) main = "Graph structure"
  	g_igraph <- graph.adjacency( true_graph, mode = "undirected", diag = FALSE )
	
    plot.igraph( g_igraph, main = main, layout = layout, ... )
}		
       
