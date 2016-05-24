

##########################################################################
# alternating least squares for Haberman linking
.linking.haberman.als <- function(logaM , wgtM , maxiter , conv,
	   progress , est.type){
	#****************
    iter <- 0
	parchange <- 1000
	NS <- ncol(logaM)
	NI <- nrow(logaM)
	logaAt <- rep(0,NS)
	#*** begin algorithm
	while( ( parchange > conv ) & (iter < maxiter) ){
		logaAt0 <- logaAt
		# calculate average item slope
		logaj <- rowSums( ( logaM - matrix( logaAt,NI,NS,byrow=TRUE) ) * wgtM , na.rm=TRUE)
		# calculate adjusted mean slope
		logaMadj <- logaM - logaj
		logaAt <- colSums( logaMadj * wgtM , na.rm=TRUE ) / colSums( wgtM , na.rm=TRUE )
		logaAt[1] <- 0
		parchange <- max( abs( c( logaAt0 - logaAt ) ) )
		if (progress){
			cat( paste0( "** " , est.type , " estimation | Iteration " , iter  , " | " , 
				"Max. parameter change = " , round( parchange , 6 ) ) , "\n")
			utils::flush.console()
				}
		iter <- iter + 1
			}
		if (progress){cat("\n")	}			
	#*** end algorithm
	res <- list( "logaAt"=logaAt , "logaj" = logaj )
	return(res)
		}
##########################################################################
