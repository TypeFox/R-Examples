
##########################################################
IRTLikelihood.ctt <- function( y , errvar , theta = NULL ){
    N <- length(y)
	if ( is.null(theta) ){
		theta <- seq(-6,6,len=21)
			}
	TP <- length(theta)			
	res <- matrix( NA , nrow=N , ncol=TP)
    errsd <- sqrt(errvar)	
	for (tt in 1:TP){
		res[ , tt] <- stats::dnorm( y , mean=theta[tt] , sd = errsd )		
		}	
	attr(res,"theta") <- matrix( theta , ncol=1 )
	attr(res,"prob.theta") <- NA
	attr(res,"G") <- 1
	class(res) <- "IRT.likelihood"
	return(res)
		}
##########################################################		