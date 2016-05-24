
###########################################
# simulate from a Dirichlet distribution
dirichlet.simul <- function( alpha ){
    # alpha is a matrix input
    #-----
    N <- nrow(alpha)
    K <- ncol(alpha)
    ygamma <- 0*alpha
    for (ii in 1:K){   # ii <- 1
        ygamma[,ii] <- stats::rgamma( n=N , shape=alpha[,ii] )
                }
    x <- ygamma / rowSums(ygamma)
    return(x)
            }
#################################################
# derivative of digamma function
digamma1 <- function(x,h=.001){
    ( base::digamma(x+h) - base::digamma(x-h) ) / (2*h)
                }
##################################################
# Maximum likelihood estimation of distribution parameters
dirichlet.mle <- function( x ,  weights=NULL , eps=10^(-5),convcrit=.00001 , maxit=1000,
		oldfac = .3 , progress=FALSE){
	
	#***
    N <- nrow(x)
    K <- ncol(x)
    # compute log pbar
	x <- ( x+eps ) / ( 1 + 2*eps )	
	x <- x / rowSums(x) 
	N <- nrow(x)
	
	if ( is.null(weights) ){
		weights <- rep(1,N) }
	weights <- N * weights / sum( weights )		
#    log.pbar <- colMeans( log( x+eps ) )
	log.pbar <- colMeans( weights * log( x ) )
    # compute inits
#    alphaprob <- colMeans( x  )
#    p2 <- mean( x[,1]^2  )	
    alphaprob <- colMeans( x * weights )
    p2 <- mean( x[,1]^2 * weights )
    xsi <- ( alphaprob[1] - p2 ) / ( p2 - ( alphaprob[1] )^2 )
    alpha <- xsi * alphaprob 
    K1 <- matrix(1,K,K)
    conv <- 1
	iter <- 1
	#******************************
    # BEGIN iterations
    while( ( conv > convcrit ) & (iter < maxit) ){
        alpha0 <- alpha
        g <- N * base::digamma( sum(alpha ) ) - N * base::digamma(alpha) + N * log.pbar
        z <- N * digamma1( sum(alpha ))
        H <- diag( -N*digamma1( alpha ) ) + z
        alpha <- alpha0 - solve(H , g )
		alpha[ alpha < 0 ] <- 10^(-10)
		alpha <- alpha0 + oldfac*( alpha - alpha0 )
        conv <- max( abs( alpha0 - alpha ) )	
		if (progress){     print( paste( iter , sum(alpha) , conv) ) }
		iter <- iter+1
		utils::flush.console()
                }
    alpha0 <- sum(alpha)
    xsi <- alpha / alpha0
    res <- list( "alpha"=alpha , "alpha0" = alpha0 , "xsi" = xsi )
    return(res)
        }
##############################################################
    				