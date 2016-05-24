

#####################################################################
# calculate log likelihood function in logistic model
.calc.ll.logistic <- function( theta , b , freq.correct , wgt ){
    eps <- 10^(-10)
	TP <- length(theta)
	I <- length(b)
	# calculate probabilities
	prob <- stats::plogis( theta + matrix( b , nrow=TP , ncol=I, byrow=TRUE) )
	# calculate log-likelihood
	ll <- freq.correct*wgt*log( prob + eps) + (1-freq.correct)*wgt*log( 1-prob+eps )	
	res <- ll
	return(res)
	}
#####################################################################
# update theta parameter in logistic distribution
.update.theta.logistic <- function( theta , b , freq.correct , wgt , 
		numdiff.parm , max.increment){
	#*****
	h <- numdiff.parm	
	# update theta parameter
	ll0 <- .calc.ll.logistic( theta , b , freq.correct , wgt )
	ll1 <- .calc.ll.logistic( theta+h , b , freq.correct , wgt )	
	ll2 <- .calc.ll.logistic( theta-h , b , freq.correct , wgt )			
	ll0 <- rowSums(ll0)
	ll1 <- rowSums(ll1)
	ll2 <- rowSums(ll2)	
	# derivative
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
    # second order derivative
    # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    # change in item difficulty
    d2[ abs(d2) < 10^(-10) ] <- 10^(-10)
    increment <- - d1 / d2
#	ci <- ceiling( abs(increment) / ( abs( max.increment) + 10^(-10) ) )
#    increment <- ifelse( abs( increment) > abs(max.increment)  , 
#                                 increment/(2*ci) , increment )	
	increment <- ifelse( abs( increment) > abs(max.increment)  , 
									 sign(increment)*max.increment , increment )
	theta <- theta + increment
    res <- list("theta" = theta , "ll" = sum(ll0) )
	return(res)
	}
##############################################################
#####################################################################
# update b parameter in logistic distribution
.update.b.logistic <- function( theta , b , freq.correct , wgt , 
		numdiff.parm , max.increment){
	#*****
	h <- numdiff.parm	
	
		# update theta parameter
		ll0 <- .calc.ll.logistic( theta , b , freq.correct , wgt )
		ll1 <- .calc.ll.logistic( theta , b+h , freq.correct , wgt )	
		ll2 <- .calc.ll.logistic( theta , b-h , freq.correct , wgt )			
		ll0 <- colSums(ll0)
		ll1 <- colSums(ll1)
		ll2 <- colSums(ll2)	
		# derivative
		d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
		# second order derivative
		# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
		d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
		# change in item difficulty
		d2[ abs(d2) < 10^(-10) ] <- 10^(-10)
		increment <- - d1 / d2
#		ci <- ceiling( abs(increment) / ( abs( max.increment) + 10^(-10) ) )
#		increment <- ifelse( abs( increment) > abs(max.increment)  , 
#									 increment/(2*ci) , increment )	
		increment <- ifelse( abs( increment) > abs(max.increment)  , 
									 sign(increment)*max.increment , increment )	

		b <- b + increment

    res <- list("b" = b , "ll" = sum(ll0) )
	return(res)
	}
##############################################################