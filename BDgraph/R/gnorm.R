# To compute Normalizing constant of G-Wishart distribution according to ATAY-KAYIS AND MASSAM (2005)
gnorm = function( adj.g, b = 3, D = diag( ncol(adj.g) ), iter = 100 )
{
	if ( b <= 2 ) stop( "In G-Wishart distribution parameter 'b' has to be more than 2" )
	if( is.null( adj.g ) ) stop( "Adjacency matrix should be determined" )

	G <- as.matrix( adj.g )
	if( sum( ( G == 1 ) * ( G == 0 ) ) != 0 ) stop( "Elements of matrix G should be zero or one" )	

	G[ lower.tri( G, diag = TRUE ) ] <- 0

	sumrowG = rowSums(G)
	sumcolG = colSums(G)
	sumG    = sum(G)
	p       = nrow(G)

	Ti      = chol( solve(D) )
	H       = Ti / t( matrix( rep( diag(Ti) ,p ), p, p ) )

	check_H = identical( H, diag(p) ) * 1

	nu  = sumrowG
	f_T = c( rep( 0, iter ) )

	result = .C( "log_exp_mc", as.integer(G), as.integer(nu), as.integer(b), as.double(H), as.integer(check_H), 
	              as.integer(iter), as.integer(p), f_T = as.double(f_T), PACKAGE = "BDgraph" )
	f_T    = c( result $ f_T )
	
	log_Ef_T = log( mean( exp( - f_T / 2 ) ) )

	c_dT = ( sumG / 2 ) * log( pi ) + ( p * b / 2 + sumG ) * log( 2 ) +
	          sum( lgamma( ( b + sumrowG ) / 2 ) ) + sum( ( b + sumrowG + sumcolG ) * log( diag(Ti) ) )
	  
	logIg = c_dT + log_Ef_T
	
	return( logIg ) 
}
      
