ini.lambda <- function(x, y, k, w, xi){
	A <- as.matrix(eigen(var(x))$vectors[ ,1:k])
	y_star <- x %*% A
	lambda_beta_max <- 0
	for( j in 1:ncol(y_star) )
	{
		for( l in 1:ncol(x) )
		{
			s <- t(x[ , l ]) %*% y_star[ , j ] * 2 * w / ( 1 - xi )
			if( lambda_beta_max < abs( s ) ) lambda_beta_max <- abs( s )
		}
	}
	ss <- lambda_beta_max - lambda_beta_max*10^(-5)
	return(ss)
}
