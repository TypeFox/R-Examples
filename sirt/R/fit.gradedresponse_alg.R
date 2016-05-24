

#####################################################################
# calculate log likelihood function in graded response model
.calc.ll.grm <- function( theta , b , b.cat , freq.categories){
    eps <- 10^(-10)
	TP <- length(theta)
	I <- length(b)
	K <- length(b.cat)	
	prob1 <- prob <- array( 1 , dim=c(TP,I,K+1) )
	for (kk in 1:K){
		prob1[,,kk+1] <- stats::plogis( theta + matrix( b , nrow=TP , ncol=I, byrow=TRUE) + b.cat[kk] )
		prob[,,kk] <- prob1[,,kk]-prob1[,,kk+1]
				 }
	kk <- K+1
	prob[,,kk] <- prob1[,,kk]
	prob[ prob < eps ] <- eps
	# calculate log-likelihood
	ll <- freq.categories * log( prob )
	res <- list("ll"=ll , "prob"=prob )	
	return(res)
	}
#####################################################################
# update theta parameter in logistic distribution
.update.theta.grm <- function( theta , b , b.cat , freq.categories ,
		numdiff.parm , max.increment){
	#*****
	h <- numdiff.parm	
	# update theta parameter
	ll0 <- .calc.ll.grm( theta , b , b.cat , freq.categories)
	prob.grm <- ll0$prob
	ll0 <- ll0$ll	
	ll1 <- .calc.ll.grm( theta+h , b , b.cat , freq.categories)$ll	
	ll2 <- .calc.ll.grm( theta-h , b , b.cat , freq.categories)$ll	
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
    res <- list("theta" = theta , "ll" = sum(ll0) , "prob.grm"=prob.grm )
	return(res)
	}
#####################################################################
# update b parameter in logistic distribution
.update.b.grm <- function( theta , b , b.cat , freq.categories ,
		numdiff.parm , max.increment){
	#*****
	h <- numdiff.parm	
	# update theta parameter
	ll0 <- .calc.ll.grm( theta , b , b.cat , freq.categories)
	ll0 <- ll0$ll	
	ll1 <- .calc.ll.grm( theta , b+h , b.cat , freq.categories)$ll	
	ll2 <- .calc.ll.grm( theta , b-h , b.cat , freq.categories)$ll				
    ll0 <- rowSums( colSums(ll0))
    ll1 <- rowSums( colSums(ll1))
    ll2 <- rowSums( colSums(ll2))	
	# derivative
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
    # second order derivative
    # f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
    d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
    # change in item difficulty
    d2[ abs(d2) < 10^(-10) ] <- 10^(-10)
    increment <- - d1 / d2
	increment <- ifelse( abs( increment) > abs(max.increment)  , 
									 sign(increment)*max.increment , increment )
	b <- b + increment
    res <- list("b" = b , "ll" = sum(ll0) )
	return(res)
	}
##############################################################
#####################################################################
# update b parameter in logistic distribution
.update.bcat.grm <- function( theta , b , b.cat , freq.categories ,
		numdiff.parm , max.increment){
	#*****
	h <- numdiff.parm		
	b.catN <- 0*b.cat
	for (kk in seq(1,length(b.cat))){
		e1 <- b.catN
		e1[kk] <- 1
		# update theta parameter
		ll0 <- .calc.ll.grm( theta , b , b.cat , freq.categories)
		ll0 <- ll0$ll	
		ll1 <- .calc.ll.grm( theta , b , b.cat+h*e1 , freq.categories)$ll	
		ll2 <- .calc.ll.grm( theta , b , b.cat-h*e1 , freq.categories)$ll				
		ll0 <- sum(ll0)
		ll1 <- sum(ll1)
		ll2 <- sum(ll2)
		# derivative
		d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
		# second order derivative
		# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
		d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2
		# change in item difficulty
		d2[ abs(d2) < 10^(-10) ] <- 10^(-10)
		increment <- - d1 / d2
		increment <- ifelse( abs( increment) > abs(max.increment)  , 
										 sign(increment)*max.increment , increment )
		b.cat[kk] <- b.cat[kk] + increment
			}
    res <- list("b.cat" = b.cat , "ll" = sum(ll0) )
	return(res)
	}
##############################################################