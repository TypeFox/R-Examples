
#####################################################
# rm.smooth.distribution
rm.smooth.distribution <- function( theta.k , pi.k , est.mean = FALSE ,
			skillspace="normal" ){
		m2 <- 0
		if ( est.mean ){
			m2 <- sum( theta.k * pi.k )
						}
		w2 <- sum( theta.k^2 * pi.k ) - m2^2
		sigma <- sqrt(w2)
		if ( skillspace == "normal" ){
			pi.k <- stats::dnorm( theta.k , mean= m2 , sd=sigma )
			pi.k <- pi.k / sum( pi.k )
							}
		res <- list( "mu"= m2 , "sigma" = sigma , "pi.k"=pi.k)
		return(res)
			}
#######################################################			