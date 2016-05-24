
#########################################################
# linear approximation of confirmatory factor analysis 
IRT.linearCFA <- function( object , group=1 ){
		# preliminaries
		expc <- IRT.expectedCounts( object )
		n.ik <- expc[,,,group]
		theta <- attr( expc , "theta" )
		pi.k <- attr( expc , "prob.theta" )
		I <- dim(n.ik)[1]
		K <- dim(n.ik)[2]
		TP <- nrow(theta)
		D <- ncol(theta)
		colnames(theta) <- paste0( "theta.Dim" , 1:D)
		# moments of trait distribution
		M.trait <- rep(0,D)
		SD.trait <- rep(0,D)
		for (dd in 1:D){
			M.trait[dd] <- m0 <- stats::weighted.mean( theta[,dd]  , pi.k )
			sd0 <- stats::weighted.mean( theta[,dd]^2  , pi.k )
			SD.trait[dd] <- sqrt( sd0 - m0^2 )
						}
		# output data frame
		dfr <- as.data.frame( matrix( 0 , nrow=I , ncol= 1+2+1+D+2 ) )
		colnames(dfr) <- c("item" , "Mlat" , "SDlat" , "mu" , 
					paste0("load.Dim",1:D) , "ResidVar" , "h2")
		dfr$item <- dimnames(expc)[[1]]

		#*************************************
		# loop over items ii		
		for (ii in 1:I){
			# ii <- 1  # Item ii
			exp.ii <- n.ik[ii,,]


			# create data frame for linear modelling
			dfr.ii <- NULL

			for (kk in 1:K){		# loop categories
				# kk <- 1
				dfr.ii.kk <- data.frame( "cat" = kk - 1, "count" = exp.ii[kk,] )
				dfr.ii.kk <- cbind( dfr.ii.kk , theta )
				dfr.ii <- rbind( dfr.ii , dfr.ii.kk )
						}

			# linear approximation factor model
			form <- paste0( "cat ~ " , paste0( colnames(theta) , collapse = " + " ) )
			mod <- stats::lm( stats::as.formula(form) , data = dfr.ii , weights = dfr.ii$count )
			# print( summary(mod) )

			# fitted values
			fitted_mod <- stats::fitted(mod)
			resid_mod <- stats::resid(mod)
			dfr[ii,"ResidVar"] <- stats::weighted.mean( resid_mod^2 , dfr.ii$count )
			dfr[ii, c("mu" , paste0("load.Dim" , 1:D) ) ] <- coef(mod)

			# latent mean and latent SD (model-implied)
			dfr[ ii , "Mlat" ] <- M1 <- stats::weighted.mean( dfr.ii$cat , dfr.ii$count )
			V1 <- stats::weighted.mean( dfr.ii$cat^2 , dfr.ii$count )
			dfr[ ii , "SDlat" ] <- sqrt( V1 - M1^2 )		
			# communality
			dfr[ ii , "h2" ] <- 1 - ( dfr[ii, "ResidVar"] ) /  dfr[ii,"SDlat"]^2 

			# standardized loadings
			dfr2 <- data.frame( "item" = dfr$item )
			
			for (dd in 1:D){
				# dd <- 1
				dfr2[ , paste0("stand.load.Dim" , dd) ] <- dfr[ , paste0("load.Dim" , dd) ] * 
										SD.trait[dd] / dfr[ , "SDlat" ]
							}
						}
		#*****************************************				
		# output				
		res <- list( "loadings" = dfr , "stand.loadings" = dfr2 , 
						"M.trait" = M.trait , "SD.trait" = SD.trait )
		class(res) <- "IRT.linearCFA"
		return(res)
				}
############################################################

				
##########################################################				
# summary method
summary.IRT.linearCFA <- function( object , ... ){
	cat("Unstandardized Loadings\n\n")
	obji <- object$loadings
	for (vv in 2:(ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
						}
	print(obji)
	cat("\nStandardized Loadings\n\n")
	obji <- object$stand.loadings
	for (vv in 2:(ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
						}
	print(obji)											
			}
##############################################################