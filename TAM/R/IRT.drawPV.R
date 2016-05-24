

###################################################
# draw plausible values
IRT.drawPV <- function( object , NPV = 5 ){
		post1 <- IRT.posterior( object )
		theta <- attr( post1 , "theta" )
		hwt <- as.matrix(post1)
		rownames(hwt) <- NULL
		attr(hwt,"theta") <- NULL
		attr(hwt,"prob.theta") <- NULL
		N <- nrow(hwt)
		TP <- nrow(theta)
		D <- ncol(theta)
		pvmatr <- matrix( NA , nrow=N , ncol=NPV*D )
		l1 <- NULL
		for (pp in 1:NPV){
			l1  <- c( l1 , paste0("PV" , pp , ".Dim" ,1:D) )
						}
		colnames(pvmatr) <- l1

		#********************************************
		# draw plausible values for uni-dimensional models
		if (D==1){
			theta <- matrix( theta , nrow=N , ncol=TP , byrow=TRUE )
		
			m1 <- rowSums( hwt * theta )
		    sd1 <- sqrt( rowSums( hwt * theta^2 ) - m1^2 )
			for (pp in 1:NPV){
				pvmatr[,pp] <- stats::rnorm( N , mean = m1 , sd = sd1 )
							}    
				}
		#*************************************************
		# draw plausible values for multidimensional models
		if (D>1){
			Sigma <- matrix( 0 , nrow=D , ncol=D)
			mu <- rep(0,D)
			
			for (nn in 1:N){
				# nn <- 1
				hwt.nn <- hwt[nn, ]
				for (dd in 1:D){
					mu[dd] <- sum( hwt.nn * theta[,dd] )
								}
				for (dd1 in 1:D){
					for (dd2 in dd1:D){			
						# dd1 <- 1
						# dd2 <- 1
						Sigma[dd1,dd2] <- sum( hwt.nn * theta[,dd1] * theta[,dd2] ) 
						Sigma[dd1,dd2] <- Sigma[dd1,dd2] - mu[dd1]*mu[dd2]
						Sigma[dd2,dd1] <- Sigma[dd1,dd2]
										}
								}
				pv_nn <- MASS::mvrnorm( NPV , mu = mu , Sigma = Sigma )				
				pvmatr[nn, ] <- as.vector( matrix( t(pv_nn) , nrow=NPV*D , ncol=1,byrow=TRUE ))
							}
				}
		#************ end imputation function						
		return(pvmatr)

				}