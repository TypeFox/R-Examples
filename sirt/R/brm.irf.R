#############################################################
# item response function (discretized) beta response model
brm.irf <- function( Theta , delta , tau , ncat , thdim=1 , eps=1E-10 ){
	TP <- nrow(Theta)
	K <- ncat
	# compute mid points
	mp <- seq( 1 / (2*(ncat) ) , 1 , 1/ncat  )
	probs <- matrix( 0 , nrow=TP , ncol=K )
	eps <- 1E-10
	# compute beta shape parameters
	m1 <- Theta[,thdim] - delta + tau
	m1 <- exp( m1 / 2 )
	m2 <- - Theta[,thdim] + delta + tau
	m2 <- exp( m2 / 2 )
	for (cc in 1:ncat){
		# cc <- 1
		probs[,cc] <- stats::dbeta( mp[cc] , shape1 = m1 , shape2 = m2 )
					   }
	probs <- probs + eps
	probs <- probs / base::rowSums(probs)
	return(probs)
		}
################################################################