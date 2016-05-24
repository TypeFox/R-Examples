###################################################
# auxiliary function reduced skill space
gdina.reduced.skillspace <- function( ntheta , Z , 
	reduced.skillspace.method=2 , eps=10^(-10) ){		
		#***********************************
		ntheta <- ntheta / sum(ntheta)
		lntheta <- matrix(log(ntheta+eps),ncol=1 )
		V <- diag( ntheta)		
		#---------------------------
		#*** skill space method 1 (CDM <= 2.5)
		if ( reduced.skillspace.method == 1){
			Z1 <- crossprod(Z , V ) %*% Z
			diag(Z1) <- diag(Z1)+eps
			covbeta <- solve( Z1 )   
			beta <- covbeta  %*% ( crossprod(Z , V ) %*% lntheta )		
							}
		#------------------------------
		#*** skill space method 2 (CDM >= 2.6)
		if ( reduced.skillspace.method == 2){		
			mod <- stats::lm( lntheta ~ 0 + Z , weights = ntheta )
			beta <- matrix( mod$coef , nrow=ncol(Z) , ncol=1 )
			beta[ is.na(beta ) ] <- 0
						}
		#*******************************************
		pred.ntheta <- exp( Z %*% beta )
		# calculate attribute probability
		attr.prob <- ( pred.ntheta / sum(pred.ntheta ) )[,1]
		#***** output
		res <- list("beta"=beta , "pred.ntheta"=pred.ntheta , 
				"attr.prob"=attr.prob)
		return(res)		
					}
############################################################					
