
###################################################
# plot item response functions
# fitted object for which the class IRT.irfprob exists
IRT.irfprobPlot <- function( object , items=NULL , 
       min.theta=-4 , max.theta=4 , cumul=FALSE , 
	   smooth=TRUE , ask=TRUE ,
	   n.theta = 40 , package = "lattice" ,  ... ){	
	#************************************   
	# extract item response functions
	irfprob <- IRT.irfprob( object )
    irfprob[ is.na(irfprob) ] <- 0
	theta <- attr( irfprob , "theta" )
	D <- ncol(theta)
	
	if ( ! is.null( items) ){
		irfprob <- irfprob[ items ,,,drop=FALSE ]
							}
	

	items.labels <- dimnames(irfprob)[[1]]
	I <- length(items.labels)
	#***********************************
	# theta grid
	theta.grid <- seq( min.theta , max.theta , length=n.theta )
	
	#***********************************
	#**** plot
	
	for (ii in 1:I){
	# ii <- 1
	
		irf.ii <- irfprob[ ii ,,]
		# compute maximum number of categories
		rs <- rowSums( irf.ii , 1 , na.rm=TRUE )
		K <- sum( rs > 0 )
		vv <- 0
		for (dd in 1:D){
			# dd <- 1
			a1 <- stats::aggregate( irf.ii[2,] , list(theta[,dd]) , mean )		
			if ( stats::sd(a1[,2])> 1E-7 ){
				vv <- dd
						}		
						}
		dd <- vv					
	
		# compute functions
		irf1 <- t( irf.ii )
		theta1 <- theta[,dd]
		a1 <- stats::aggregate( irf1 , list( theta1 ) , mean )
		theta1 <- a1[,1]	
		dfr <- NULL
		dfr1 <- data.frame( "theta" = theta.grid )
		# btheta <- bs( theta.grid )
		# create data frame for plot in lattice
		hh <- 1
		if ( cumul ){ hh <- 2 }
		
		for (kk in hh:K){
			# kk <- 1
			v1 <- as.numeric(a1[,kk+1])
			if ( cumul ){
			   v1 <- rowSums( a1[ , seq( kk+1 , K+1) , drop=FALSE] )
						}
			if ( smooth ){			
				eps <- 1E-5
				v1 <- ( v1 + eps ) / ( 1 + 2*eps )				
				v1 <- stats::qlogis(v1)
				dat <- data.frame( "theta" = theta1 , "y" = v1 )		
				mod <- stats::loess( y ~ theta  , data=dat)
				ypred <- stats::plogis( stats::predict( mod , dfr1 ) )
						} else {
				ypred <- v1
				theta.grid <- theta1
						}
			dfr1 <- data.frame( "theta" = theta.grid , "cat" = kk - 1 , "prob" = ypred )
			dfr <- rbind( dfr , dfr1 )
						}
									
						
		item <- items.labels[ii]
		main <- paste0('Trace lines for item ', item)
		vkey <- paste0("Cat " , 0:(K-1) )
		L1 <- K
		#**************************************
		# package lattice
		if ( package == "lattice"){		
			print( 
					lattice::xyplot( prob ~ theta , data=dfr , group=cat , ylim=c(-.1 , 1.1) , type="o" ,  
						auto.key= TRUE  ,       
						par.settings = list(superpose.symbol = list(pch = 1:L1))	,			
						ylab = expression(P(theta)), xlab = expression(theta) , 
						main=main , lty=1:L1 , pch=1:L1 ,
						...
								)  
							)
					}
		#******************************************
		# package graphics
		if ( package == "graphics" ){
		       kk <- 0
		       dfr1a <- dfr[ dfr$cat == kk , ] 
				graphics::plot( dfr1a$theta , dfr1a$prob , ylim=c(-.1,1.1) , 	
						ylab = expression(P(theta)), xlab = expression(theta) ,
						col=kk+2 , pch = kk+1 , type="o" , main=main , ...
									)
		      for (kk in seq(1,K-1) ){
				dfr1a <- dfr[ dfr$cat == kk , ]		
				graphics::lines( dfr1a$theta , dfr1a$prob , pch=kk+1 , col=kk+2 )					
				graphics::points( dfr1a$theta , dfr1a$prob , pch=kk+1 , col=kk+2 )	
									
									}
			graphics::legend( min( dfr1a$theta) , 1.1 , vkey , pch = 1:K , col= 1 + 1:K ,
						horiz =TRUE , lty= 1)
									}
					
		graphics::par( ask=ask )	

			}

		}
######################################################			