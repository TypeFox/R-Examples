


#######################################################################
# moments of distribution
.mml.3pl.distributionmoments <- function( D , G , pi.k , theta.k ){

	D <- ncol(theta.k)
	if ( is.vector( theta.k) ){  
			D <- 1 
			theta.k <- matrix( theta.k , ncol=1 )
				}

	mean.trait <- sd.trait <- skewness.trait <- matrix( 0 , nrow=D , ncol=G )
	for (dd in 1:D){
		for (gg in 1:G){
		mean.trait[dd,gg] <- sum( theta.k[,dd] * pi.k[ , gg ] )
		sd.trait[dd,gg] <- sqrt( sum( theta.k[,dd]^2 * pi.k[ , gg ] ) - mean.trait[dd,gg]^2 )
		skewness.trait[dd,gg] <- sum( ( theta.k[,dd] - mean.trait[dd,gg] )^3 * pi.k[ , gg ] ) /
						sd.trait[dd,gg]^3
					}
			}
	rownames(skewness.trait) <- rownames(sd.trait) <- 
				rownames(mean.trait) <- colnames(theta.k)
	colnames(skewness.trait) <- colnames(sd.trait) <- 
				colnames(mean.trait) <- paste0("Group",1:G)
	#*****
	# correlation matrices
	correlation.trait <- as.list(1:G)
	names(correlation.trait) <- colnames(mean.trait)
	for (gg in 1:G){
		# gg <- 1
		mean.gg <- rep(0,D)
		Sigma.gg <- diag(0,D)	
		for (dd in 1:D){
			mean.gg[dd] <- sum( pi.k[,gg] * theta.k[,dd] )
				}
		for (dd1 in 1:D){
			for (dd2 in dd1:D){
		#		dd1 <- 1 ; 	dd2 <- 1
				Sigma.gg[dd1,dd2] <- sum( pi.k[,gg] * (theta.k[,dd1] - mean.gg[dd1] )*(theta.k[,dd2] - mean.gg[dd2] ) ) 
#				Sigma.gg[dd1,dd2] <- Sigma.gg[dd1,dd2] - mean.gg[dd1] * mean.gg[dd2]
				Sigma.gg[dd2,dd1] <- Sigma.gg[dd1,dd2]
									}
						}
		rownames(Sigma.gg) <- colnames(Sigma.gg) <- rownames(mean.trait)
		correlation.trait[[gg]] <- stats::cov2cor(Sigma.gg + diag(10^(-5),D) )
					}	
	res <- list( "mean.trait"=mean.trait , "sd.trait" = sd.trait , 
				"skewness.trait" = skewness.trait , "correlation.trait"=correlation.trait)
    return(res)				
		}