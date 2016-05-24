

#############################################################################
# Robust linking
linking.robust <- function(  itempars ){
	#*******************************
	itempars0 <- itempars
    itempars <- na.omit(itempars)
    pars <- itempars[,2] - itempars[,3]
    items <- paste(itempars[,1])    
	names(pars) <- items
    I <- length(pars)

	x <- sort(pars)
	kvec <- seq(1 , floor( (I-1)/2  ) )
	KK <- length(kvec)
	se <- meanpars <- rep(NA, KK )
	# define trimming factor
	for (kk in 1:KK){
		# arrange calculations
		N <- length(x)
		k <- kk
		indkk <- seq( k+1 ,  N-k ,1)
		x0 <- x[ indkk ]
		# compute winsorized mean
		trim.mean <- mean( x0 )
		swk2 <- k * ( x[ k] - trim.mean ) ^2 + sum( ( x0 - trim.mean )^2 ) + k * ( x[ N - k + 1] - trim.mean )^2
		# standard error
		se.trimmean <- sqrt( swk2 ) / sqrt( (N-2*k) * ( N - 2*k - 1 ) )
		# output
		meanpars[kk] <- trim.mean
		se[kk] <- se.trimmean
				}
            
		v1 <- paste0("k" , 0:KK)
		meanpars <- c( mean(x) , meanpars )
		se <- c( stats::sd(x) / sqrt(I) , se )
		names(meanpars) <- v1
		names(se) <- v1

	# arrange output
    res1 <- list()
    res1$ind.kopt <- ind.kopt <- which.min( se )
    res1$kopt <- kvec[ ind.kopt ] - 1
    res1$meanpars.kopt <- meanpars[ ind.kopt ]
    res1$se.kopt <- se[ ind.kopt ]
    res1$meanpars <- meanpars
    res1$se <- se    
    res1$sd <- stats::sd(x)
	res1$mad <- stats::mad(x)
    res1$k.robust <- c(0,kvec  )
	res1$I <- I
	res1$itempars <- itempars0
	class(res1) <- "linking.robust"
    return(res1)
        }
#############################################################################
# S3 plot method
plot.linking.robust <- function( x ,  ... ){
        graphics::par( mfrow=c(2,1))
	KK <- length(x$k.robust)
    graphics::plot( x$k.robust , x$meanpars[1:KK] , type="l" , xlab="k" , 
		ylab="Linking constant" ,
        main="Linking constant")
    graphics::points( 0 , x$meanpars[1] , pch=16 , col=3 , cex=1.4 )              	
	graphics::points( x$kopt , x$meanpars.kopt , pch=17 , col=2 , cex=1.4 )              	
	#****
    graphics::plot( x$k.robust , x$se[1:KK] , type="l" , 
			main= paste0( "Standard error of linking constant (k_opt = " , round(x$kopt , 3 ),")" ) ,
            xlab="k" , ylab="Standard error")
	graphics::points( 0 , x$se[1] , pch=16 , col=3 , cex=1.4 ) 
    graphics::points( x$kopt , x$se.kopt , pch=17 , col=2 , cex=1.4 )  	
        graphics::par( mfrow=c(1,1) )          
        }
#################################################################################
# S3 summary method
summary.linking.robust <- function( object , ... ){
    kmax <- length(object$k.robust)
    cat("Robust linking with trimmed mean\n\n")
    cat( paste0("Number of items = " , object$I , "\n" ) )
#    cat( paste0("Number of jackknife units = " , object$JJ , "\n\n" ) )
    cat(paste0( "Optimal trimming parameter k = " , round( object$kopt , 4 ) , " | "))	
	cat(paste0( " non-robust parameter k = " , 0 , " \n"))
    cat(paste0( "Linking constant = " , round( object$meanpars.kopt , 4 ) , " | "))
    cat(paste0( " non-robust estimate = " , round( object$meanpars[ 1 ] , 4 ) , " \n"))
    cat(paste0( "Standard error = " , round( object$se.kopt , 4 ) , " | "))
    cat(paste0( " non-robust estimate = " , round( object$se[1] , 4 ) , " \n"))
    cat(paste0( "DIF SD: MAD = " , round( object$mad , 4 ) , " (robust) | "))
    cat(paste0( "SD = " , round( object$sd , 4 ) , " (non-robust) \n"))	
#    cat(paste0( " non-robust estimate = " , round( object$sumweight[ kmax ] , 4 ) , " \n"))	
        }