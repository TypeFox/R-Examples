


#########################################################################################				 
.mstep.mml.npirt <- function( pjk , r.jk , n.jk , theta.k , npformula , npmodel, G , I ,
			npirt.monotone , ICC_model_matrix ){				 
					rjk0 <- r.jk[,,1]
					njk0 <- n.jk[,,1]
					if (G > 1){
						for (gg in 2:G){
							rjk0 <- rjk0 + r.jk[,,gg]
							njk0 <- njk0 + n.jk[,,gg]						
									}
								}
									
		if ( is.null( npformula ) ){ 
					pjk <- t(rjk0 / njk0 )
					if (npirt.monotone){
					    # monotone smoothing
						pjk <- monoreg.colwise(yM=pjk, wM=t(njk0) )
									}
					}

		# estimation using a formula for ICC estimation
		if ( ! is.null(npformula) ) {
				LK <- length(theta.k)
						cat("ICC estimation |")	
						prbar <- floor( 10 * ( 1:I )	/ I )
						prbar <- c( 1 , diff(prbar))
						for (ii in 1:I){		
								#ii <- 3
#							dfr1 <- data.frame( "theta" = theta.k , "y" = 1 , "wgt" = rjk0[ii,] )
#							dfr0 <- data.frame( "theta" = theta.k , "y" = 0 , "wgt" = njk0[ii,] - rjk0[ii,] )
#							dafr <- data.frame( rbind( dfr0 , dfr1 ) )
							if ( prbar[ii] == 1){ cat("~"); utils::flush.console() }
#							theta <- dafr$theta.k
#							wgt <- dafr$wgt
#							y <- dafr$y
#							ICC_ <- model.matrix( npformula[[ii]] , dafr )						
							y <- rep( c(0,1) , each=LK)
							wgt <- c( njk0[ii,] - rjk0[ii,] , rjk0[ii,]  )
							ICC_ <- ICC_model_matrix[[ii]]
							npmodel[[ii]] <- stats::glm( y ~ 0 + ICC_ , weights = wgt , family="binomial" ,
										control=list(maxit=4) )						
#							npmodel[[ii]] <- glm( npformula[[ii]] , 
#										data = dafr , weights = dafr$wgt , family="binomial")					
							pjk[,ii] <- stats::fitted( npmodel[[ii]] )[ seq( 1 , LK ) + LK ]
								}
							cat("\n")		
						}
					res <- list( "pjk" = pjk , "npmodel"=npmodel)				
							}
##########################################################################################

