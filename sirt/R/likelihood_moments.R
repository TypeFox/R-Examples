
##############################################################
likelihood_moments <- function( likelihood , theta=NULL ){		
		if ( is.null(theta) ){
			theta <- attr( likelihood , "theta" )
							}		
		nstud <- nrow(likelihood)
		TP <- ncol(likelihood)
		thetaM <- matrix( theta , nstud , TP , byrow=TRUE )
#		w1 <- rowSums( likelihood )
#		M1 <- rowSums( thetaM * likelihood ) / w1	
#		SD1 <- rowSums( thetaM^2 * likelihood ) / w1

		likelihood <- likelihood / rowSums( likelihood )
		M1 <- rowSums( thetaM * likelihood )
		SD1 <- rowSums( thetaM^2 * likelihood )
		SD1 <- sqrt( SD1 - M1^2 )	
		res <- list( M = M1 , "SD" = SD1 )
		return(res)
				}
##############################################################				
		
	