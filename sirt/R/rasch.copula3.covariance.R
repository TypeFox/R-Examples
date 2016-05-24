

######################################
# estimation of covariance
rasch.copula3.covariance <- function( f.qk.yi , Sigma , theta.k , N ,
		mu.fixed , variance.fixed , D , est.corr , irtmodel ,
		freqwgt ){
		Sigma.cov <- Sigma
		delta.theta <- 1
		indD <- rep( 1:nrow(f.qk.yi) , freqwgt )
		hwt <- f.qk.yi[ indD , ]	
		N <- sum( freqwgt)
#		if (qmcnodes){ 
#			hwt <- hwt / nrow(theta.k) 
#			hwt <- hwt / rowSums( hwt )
#				}	
		thetabar <- hwt%*%theta.k
		mu <- rep(0,D)		
		# calculation of mu
#		mu <- colSums( thetabar ) / N
#		if ( ! is.null(mu.fixed ) ){
#			  if (is.matrix(mu.fixed) ){	 
#			    mu0 <- mu
#				mu[ mu.fixed[,1] ] <- mu.fixed[,2]
#									}
#							}
		# calculation of the covariance matrix
#		theta.k.adj <- theta.k - matrix( mu , nrow=nrow(theta.k) , 
#									ncol=ncol(theta.k) , byrow=TRUE)							
		theta.k.adj <- theta.k 
		for (dd1 in 1:D){
			for (dd2 in dd1:D){
				tk <- theta.k.adj[,dd1]*theta.k.adj[,dd2]
				h1 <- ( hwt %*% tk ) * delta.theta			
				Sigma.cov[dd1,dd2] <- sum( h1  ) / N
				if (dd1 < dd2 ){ Sigma.cov[dd2,dd1] <- Sigma.cov[dd1,dd2] }
										
									}
								}
		Sigma.cov0 <- Sigma.cov
								
		if ( est.corr ){ Sigma.cov <- stats::cov2cor(Sigma.cov ) }					
		if ( ! is.null(variance.fixed ) ){
				Sigma.cov[ variance.fixed[,1:2,drop=FALSE] ] <- variance.fixed[,3]
				Sigma.cov[ variance.fixed[,c(2,1),drop=FALSE] ] <- variance.fixed[,3]		
#										
									}
		diag(Sigma.cov) <- diag(Sigma.cov) + 10^(-10)

		
		pi.k <- matrix( mvtnorm::dmvnorm( theta.k , mean = mu , sigma = Sigma.cov )	, ncol=1 )		
		pi.k <- pi.k / sum( pi.k )			
		res <- list( "mu"=mu , "Sigma"=Sigma.cov , "pi.k"= pi.k ,
			"Sigma0" = Sigma.cov0)				
		return(res)
					}