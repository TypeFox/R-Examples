
###########################################################################
# IRT.truescore
IRT.truescore <- function( object , iIndex = NULL , theta=NULL, Q = NULL){
		irf <- IRT.irfprob(object)
		theta0 <- attr( irf , "theta" )[,1]
		dim_irf <- dim(irf)
		if ( is.null(iIndex) ){
			iIndex <- seq( 1 , dim_irf[1] )
						}
		irf <- irf[ iIndex ,, ]
		K <- dim_irf[2] - 1
		irf[ is.na(irf)] <- 0
		I <- dim(irf)[1]
		
		if ( is.null(Q) ){
			Q <- matrix( seq(0,K) , nrow=I , ncol=K+1 , byrow=TRUE)
				}
		#*******
		# compute expected scores
		irf_K <- 0+0*irf
#		for (kk in 1:K){
#			irf_K[,kk+1,] <- kk
#						}
		for (ii in 1:I){
			irf_K[ii,,] <- as.vector(Q[ii,])
				}
				
				
		TP <- dim_irf[3]
		I <- dim(irf)[1]	
		vec <- rep(0,TP)

		for (ii in 1:I){
			# ii <- 1
			vec <- vec +  colSums( irf[ii,,] * irf_K[ii,,]  )
						}
		dfr <- data.frame( "theta" = theta0 , "truescore" = vec)
		#*******
		if ( ! is.null(theta) ){
		    ind <- sum( ! ( theta >= min(theta0) ) & ( theta <= max(theta0) ) ) > 0
			if (ind){
				h1 <- "The true score cannot be computed for some theta values.\n"
				warning(paste0( h1 , "  Use a larger theta grid in the fitted model.\n") )
					}
			theta <- theta[ ( theta >= min(theta0) ) & ( theta <= max(theta0) ) ]
			TP1 <- length(theta)
			v1 <- 1
			v2 <- TP
			for (tt in 1:TP){
				v1 <- ifelse( theta > theta0[tt] , tt , v1 )	
				v2 <- ifelse( theta < theta0[TP-tt+1] , TP-tt+1 , v2 )
							}
			vec <- dfr$truescore[ v1 ] +  ( dfr$truescore[ v2 ] - dfr$truescore[ v1 ] ) *
						( theta - dfr$theta[ v1 ] ) / ( dfr$theta[ v2 ] - dfr$theta[ v1 ] )
			dfr <- data.frame( "theta" = theta , "truescore" = vec)						
							}
		return(dfr)						
			}
###########################################################################			