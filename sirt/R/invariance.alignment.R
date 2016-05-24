################################################################
invariance.alignment__sirt045 <- function( lambda , nu , wgt ,
   align.scale=c(1,1) , align.pow=c(1,1) , eps=.0001 , 
   h= .001 ,max.increment=.2 , increment.factor=1.001 ,  maxiter =3000 ,
   conv=.0001 , fac.oldpar = 0.85  , psi0.init=NULL , alpha0.init=NULL ,
   progress=TRUE ){
    s1 <- Sys.time()
	max.incr0 <- max.increment
	#****
	G <- nrow(lambda)   # number of groups
	I <- ncol(lambda)   # number of items
	wgt <- G * wgt / colSums(wgt) 

	# missing indicator matrix: 1 - no missings
	missM <- (1-is.na(lambda))+ (1- is.na(nu))
	wgt <- wgt * missM
	lambda[ missM == 0 ] <- mean( lambda , na.rm=TRUE )
	nu[ missM == 0 ] <- mean( nu , na.rm=TRUE )
	
	group.combis <- t( combn( G , 2 ) )
	group.combis <- rbind( group.combis , group.combis[,c(2,1) ] )
	group.combis <- group.combis[ order( group.combis[,1] ) , ]

	# initial lambda0 and nu estimates for all items
	# lambda0 <- colSums( lambda * wgt ) / colSums( wgt )
	# nu0 <- colSums( nu * wgt ) / colSums( wgt )

	# initial estimates means and SDs
#	alpha0 <- rep( 0 , G )
#	psi0 <- rep(1,G)
	psi0 <- rowMeans( lambda )
	psi0 <- psi0 * ( prod( psi0 ) )^( -1/G )
	alpha0 <- rowMeans( nu )
	alpha0 <- alpha0 - alpha0[1]
	if ( ! is.null( psi0.init) ){ psi0 <- psi0.init }
	if ( ! is.null( alpha0.init) ){ alpha0 <- alpha0.init / psi0 }
	
	fopt.history <- matrix(NA , nrow=maxiter , ncol=2)
	psi0b <- psi0
	iter <- 1
	minval <- fopt <- 10^300
	fopt_change <- 10000
	Niter <- rep(NA , 2 )
	miniter0 <- fopt0 <- minval0 <- rep(NA,2)
	#################################
	# begin iterations (lambda)
	while( ( iter <= maxiter ) & ( fopt_change > conv) ){ 
		alpha0_old <- alpha0
		psi0_old <- psi0
		fopt_old <- fopt
		#*****
		# optimization group SDs        
		psi0b <- psi0
		# f(x)
		flambda <- ll0a <- align.optim.lambda( lambda=lambda , psi0=psi0 , psi0b=psi0b ,
					align.scale=align.scale[1] , align.pow=align.pow[1] ,wgt, eps=eps,group.combis)
		fopt <- ll0 <- ll0a
		
		# f(x+h)
		ll1 <- align.optim.lambda( lambda=lambda , psi0=psi0+h , psi0b=psi0b ,
					align.scale=align.scale[1] , align.pow=align.pow[1] , wgt,eps=eps,group.combis)
		
		# f(x-h)
		ll2 <- align.optim.lambda( lambda=lambda , psi0=psi0-h , psi0b=psi0b ,
					align.scale=align.scale[1] , align.pow=align.pow[1] , wgt , eps=eps,group.combis)
									
					
		# first and second derivative
		increment <- align.newton.raphson( ll0 , ll1 , ll2 , max.increment , h )
		psi0 <- psi0 + increment
		psi0[ psi0 < eps ] <- eps
		psi0 <- psi0 * ( prod( psi0 ) )^( -1/G )

		
#		fac.oldpar <- 0.5
# fac.oldpar <- 0		
		psi0 <- (1-fac.oldpar)*psi0+fac.oldpar*psi0_old
#		psi00 <- psi0
#		psi0 <- psi0_old
		#*****
		# optimization group means
#		alpha0b <- alpha0 
#		fnu <- ll0 <- align.optim.nu( lambda , nu , psi0 , psi0b ,alpha0 , alpha0b , align.scale[2] , align.pow[2] ,
#			 wgt , eps, group.combis)
#		ll1 <- align.optim.nu( lambda , nu , psi0 , psi0b ,alpha0+h , alpha0b , align.scale[2] , align.pow[2] ,
#			 wgt , eps , group.combis)
#		ll2 <- align.optim.nu( lambda , nu , psi0 , psi0b , alpha0-h , alpha0b , align.scale[2] , align.pow[2] ,
#			 wgt , eps , group.combis)
		# first and second derivative
#		increment <- align.newton.raphson( ll0 , ll1 , ll2 , max.increment , h )
#		alpha0 <- alpha0 + increment
#		alpha0 <- alpha0 - mean( alpha0 )
#		alpha0 <- alpha0 - alpha0[1]
#		psi0 <- psi00
		
#		fac.oldpar <- 0
#		alpha0 <- (1-fac.oldpar)*alpha0+fac.oldpar*alpha0_old
		
		#****
		# optimization function
#		fopt.history[iter] <- fopt <- sum( flambda + fnu )
		fopt.history[iter,1] <- fopt <- sum( flambda )
#		fopt <- sum( fopt )
		fopt_change <- abs( fopt - fopt_old )

#		alpha_change <- max( abs( alpha0 - alpha0_old ))
		psi_change <- max( abs( psi0 - psi0_old ))
		if (progress){
			cat( "---------------- ITERATION" , iter , " (LAMBDA) -----------------\n") 
			cat( "Optimization Function (lambda) =" , round( fopt ,6 ) )
			if ( iter > 1){ 
				cat( " | Absolute Change =" , round( fopt_change , 6 )) 
					} 
			cat("\n")
#			cat("** Maximum alpha parameter change =" , round( alpha_change , 6) , "\n" )
			cat("** Maximum psi parameter change   =" , round( psi_change , 6) , "\n" )
			flush.console()
					}
		if ( fopt < minval ){
			minval <- fopt
			miniter <- iter
#			alpha0.min <- alpha0
			psi0.min <- psi0
					}

		iter <- iter + 1
		Niter[1] <- iter
		max.increment <- max.increment / increment.factor
			} # end iterations
	fopt0[1] <- fopt ; minval0[1] <- minval
	miniter0[1] <- miniter
	#***********************************************************
	
	
	iter <- 1
	minval <- fopt <- 10^300
	fopt_change <- 10000
	max.incr0 -> max.increment
	psi0 <- psi0.min
	
	#################################
	# begin iterations (nu)
	while( ( iter <= maxiter ) & ( fopt_change > conv) ){ 
		alpha0_old <- alpha0
#		psi0_old <- psi0
		fopt_old <- fopt
		#*****
		# optimization group means
		alpha0b <- alpha0 
		fnu <- ll0 <- align.optim.nu( lambda , nu , psi0 , psi0b ,alpha0 , alpha0b , align.scale[2] , align.pow[2] ,
			 wgt , eps, group.combis)
		ll1 <- align.optim.nu( lambda , nu , psi0 , psi0b ,alpha0+h , alpha0b , align.scale[2] , align.pow[2] ,
			 wgt , eps , group.combis)
		ll2 <- align.optim.nu( lambda , nu , psi0 , psi0b , alpha0-h , alpha0b , align.scale[2] , align.pow[2] ,
			 wgt , eps , group.combis)
	# first and second derivative
		increment <- align.newton.raphson( ll0 , ll1 , ll2 , max.increment , h )
		alpha0 <- alpha0 + increment
		alpha0 <- alpha0 - mean( alpha0 )
		alpha0 <- alpha0 - alpha0[1]
#		psi0 <- psi00
		
#		fac.oldpar <- 0
		alpha0 <- (1-fac.oldpar)*alpha0+fac.oldpar*alpha0_old
		
		#****
		# optimization function
		fopt.history[iter,2] <- fopt <- sum( fnu )
#		fopt.history[iter] <- fopt <- sum( flambda )
#		fopt <- sum( fopt )
		fopt_change <- abs( fopt - fopt_old )

		alpha_change <- max( abs( alpha0 - alpha0_old ))
#		psi_change <- max( abs( psi0 - psi0_old ))
		if (progress){
			cat( "---------------- ITERATION" , iter , " (NU) -----------------\n") 
			cat( "Optimization Function (nu) =" , round( fopt ,6 ) )
			if ( iter > 1){ 
				cat( " | Absolute Change =" , round( fopt_change , 6 )) 
					} 
			cat("\n")
			cat("** Maximum alpha parameter change =" , round( alpha_change , 6) , "\n" )
#			cat("** Maximum psi parameter change   =" , round( psi_change , 6) , "\n" )
			flush.console()
					}
		if ( fopt < minval ){
			minval <- fopt
			miniter <- iter
			alpha0.min <- alpha0
#			psi0.min <- psi0
					}

		iter <- iter + 1
		Niter[2] <- iter
		max.increment <- max.increment / increment.factor
			} # end iterations
	fopt0[2] <- fopt ; minval0[2] <- minval ; 
	miniter0[2] <- miniter
	#***********************************************************
		
	
	
	#*****************************
	# calculate item statistics and R-squared measures
	# groupwise aligned loading
	lambda.aligned <- lambda / psi0.min
	nu.aligned <- nu - alpha0.min * lambda 
	# average aligned parameter
	itempars.aligned <- data.frame("M.lambda" = colMeans(lambda.aligned) ,
			"SD.lambda" = apply( lambda.aligned , 2 , sd , na.rm=TRUE ) ,
			"M.nu" = colMeans( nu.aligned ) ,
			"SD.nu" = apply( nu.aligned , 2 , sd , na.rm=TRUE ) 	
				)
	rownames(itempars.aligned) <- colnames(lambda)
	lambda.resid <- lambda.aligned - 
		matrix( itempars.aligned$M.lambda , G , I , byrow=TRUE )
	nu.resid <- nu.aligned - 	
		matrix( itempars.aligned$M.nu , G , I , byrow=TRUE )
		
	# R-squared measures
	Rsquared.invariance <- c(NA,NA)
	names(Rsquared.invariance) <- c("loadings" , "intercepts" )	
	expl <- psi0.min * matrix( itempars.aligned$M.lambda , G , I , byrow=TRUE)
	Rsquared.invariance["loadings"] <- 
			1 - sum( (lambda - expl)^2 , na.rm=TRUE ) / 
			sum( ( lambda - rowMeans(lambda) )^2 , na.rm=TRUE)
	expl <- matrix( itempars.aligned$M.nu , G , I , 
				byrow=TRUE) + 
			alpha0.min*matrix( itempars.aligned$M.lambda , G , I , 
				byrow=TRUE)
	Rsquared.invariance["intercepts"] <- 
			1 - sum( ( nu - expl)^2 , na.rm=TRUE ) / 
			sum( (nu - rowMeans(nu) )^2 , na.rm=TRUE)

    # correlations aligned parameters			
	rbar <- c( ai.calc.corr(t(lambda.aligned)) , ai.calc.corr(t(nu.aligned)) )			
	es.invariance <- rbind( Rsquared.invariance ,
			sqrt( 1- Rsquared.invariance ) , rbar)
	rownames(es.invariance) <- c("R2" , "sqrtU2" , "rbar")
	
	psi0 <- psi0.min
	alpha0 <- alpha0.min * psi0.min
	# original in Muthen paper: alpha / psi
	# but in this optimization alpha* = alpha / psi
	# and therefore alpha = alpha* x psi
	
	s2 <- Sys.time()
	#*****************************
	# OUTPUT:
	pars <- data.frame("alpha0"=alpha0.min*psi0.min , "psi0"=psi0.min)
	res <- list( "pars"=pars , "itempars.aligned" = itempars.aligned ,
			"es.invariance" = es.invariance , 
			"lambda.aligned" = lambda.aligned ,
			"lambda.resid" = lambda.resid ,
			"nu.aligned" = nu.aligned ,
			"nu.resid" = nu.resid ,
			"Niter" = Niter -1 , "miniter"=miniter0 ,
			"fopt"=minval0 ,
			"fopt.history" = fopt.history[1:(max(Niter)-1) , ] ,
			"align.scale"=align.scale , "align.pow"=align.pow ,
			"s1"=s1 , "s2"=s2)
	class(res) <- "invariance.alignment"
	return(res)
		}
#####################################
