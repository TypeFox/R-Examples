
##############################################################################
# Algorithm for LAMBDA
##############################################################################
inv.alignment2.lambda.alg <- function( lambda , psi0 , nu , align.scale ,
		align.pow , iter , conv , fopt_change , wgt , eps ,
		group_combis , fopt , h , max.increment , fac.oldpar , maxiter ,
		alpha0 , group.combis , G , increment.factor ){
		
	minval <- fopt <- 1E+100
	fopt_change <- 10000

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
#		flambda <- ll0a <- align.optim.lambda( lambda=lambda , psi0=psi0 , psi0b=psi0b ,
#					align.scale=align.scale[1] , align.pow=align.pow[1] ,wgt, eps=eps,group.combis)
		flambda <- ll0a <- ia_optim_lambda_R( lambda , psi0 , psi0b , align.scale[1] , align.pow[1] ,
					wgt , eps , group.combis )
				
		fopt <- ll0 <- ll0a		
		# f(x+h)
#		ll1 <- align.optim.lambda( lambda=lambda , psi0=psi0+h , psi0b=psi0b ,
#					align.scale=align.scale[1] , align.pow=align.pow[1] , wgt,eps=eps,group.combis)
		ll1 <- ia_optim_lambda_R( lambda , psi0+h , psi0b , align.scale[1] , align.pow[1] ,
					wgt , eps , group.combis )	
		# f(x-h)
#		ll2 <- align.optim.lambda( lambda=lambda , psi0=psi0-h , psi0b=psi0b ,
#					align.scale=align.scale[1] , align.pow=align.pow[1] , wgt , eps=eps,group.combis)
		ll2 <- ia_optim_lambda_R( lambda , psi0-h , psi0b , align.scale[1] , align.pow[1] ,
					wgt , eps , group.combis )														
		# first and second derivative
		increment <- align.newton.raphson( ll0 , ll1 , ll2 , max.increment , h )
		psi0 <- psi0 + increment
		psi0[ psi0 < eps ] <- eps
		psi0 <- psi0 * ( prod( psi0 ) )^( -1/G )
	
		psi0 <- (1-fac.oldpar)*psi0+fac.oldpar*psi0_old

		#****
		# optimization function
#		fopt.history[iter] <- fopt <- sum( flambda + fnu )
		fopt <- opt <- sum( flambda )
#		fopt <- sum( fopt )
		fopt_change <- abs( fopt - fopt_old )
#		alpha_change <- max( abs( alpha0 - alpha0_old ))
		psi_change <- max( abs( psi0 - psi0_old ))
		fopt <- opt
	
		if ( fopt < minval ){
			minval <- fopt
			miniter <- iter
			psi0.min <- psi0
					}

		iter <- iter + 1	
		max.increment <- max.increment / increment.factor
			} # end iterations
				
	#***********************************************************
	res <- list( "psi0" = psi0.min , "fopt" = minval )
	return(res)	
		}
#############################################################################
# auxiliary function which Calls Rcpp
ia_optim_lambda_R <- function( lambda , psi0 , psi0b , align.scale , align.pow ,
					wgt , eps , group.combis ){
     .Call("ia_optim_lambda" ,  lambda , psi0 , psi0b , align.scale , align.pow ,
					wgt , eps , group.combis-1 , PACKAGE="sirt")
											}
################################################################################

##############################################################################
# Algorithm for NU
##############################################################################
inv.alignment2.nu.alg <- function( lambda , psi0 , nu , align.scale ,
		align.pow , iter ,  conv , fopt_change , wgt , eps ,
		group_combis , fopt , h , max.increment , fac.oldpar , maxiter ,
		alpha0 , group.combis , G , increment.factor ){		
	
	minval <- fopt <- 1E+100
	fopt_change <- 10000
    psi0b <- psi0	
	
	#################################
	# begin iterations (nu)
	while( ( iter <= maxiter ) & ( fopt_change > conv) ){ 
		alpha0_old <- alpha0
		fopt_old <- fopt
		#*****
		# optimization group means
		alpha0b <- alpha0 
#		fnu <- ll0 <- align.optim.nu( lambda , nu , psi0 , psi0b ,alpha0 , alpha0b , align.scale[2] , align.pow[2] ,
#			 wgt , eps, group.combis)

		fnu <- ll0 <- ia_optim_nu_R( lambda , nu , psi0 , psi0b , alpha0 , alpha0b , align.scale[2] , align.pow[2] ,
					wgt , eps , group.combis )

					
		# ll1 <- align.optim.nu( lambda , nu , psi0 , psi0b ,alpha0+h , alpha0b , align.scale[2] , align.pow[2] ,
		# 	 wgt , eps , group.combis)
		ll1 <- ia_optim_nu_R( lambda , nu , psi0 , psi0b , alpha0+h , alpha0b , align.scale[2] , align.pow[2] ,
					wgt , eps , group.combis )			 

				
#		ll2 <- align.optim.nu( lambda , nu , psi0 , psi0b , alpha0-h , alpha0b , align.scale[2] , align.pow[2] ,
#			 wgt , eps , group.combis)
		ll2 <- ia_optim_nu_R( lambda , nu , psi0 , psi0b , alpha0-h , alpha0b , align.scale[2] , align.pow[2] ,
					wgt , eps , group.combis )				 
			 
	# first and second derivative
		increment <- align.newton.raphson( ll0 , ll1 , ll2 , max.increment , h )
		alpha0 <- alpha0 + increment
		alpha0 <- alpha0 - mean( alpha0 )
#		alpha0 <- alpha0 - alpha0[1]		
		alpha0 <- (1-fac.oldpar)*alpha0+fac.oldpar*alpha0_old
		alpha0 <- alpha0 - mean( alpha0 )
		#****
		# optimization function
		fopt <- sum( fnu )
		fopt_change <- abs( fopt - fopt_old )

		alpha_change <- max( abs( alpha0 - alpha0_old ))
		
		if ( fopt < minval ){
			minval <- fopt
			miniter <- iter
			alpha0.min <- alpha0
					}

		iter <- iter + 1
		max.increment <- max.increment / increment.factor
			} # end iterations
	
	#***********************************************************
	res <- list( "alpha0" = alpha0.min , "fopt" = minval )
	return(res)	
		}
###########################################################################		
# auxiliary function for calling Rcpp
ia_optim_nu_R <- function( lambda , nu , psi0 , psi0b , alpha0 , alpha0b , align.scale , align.pow ,
					wgt , eps , group.combis ){
	.Call( "ia_optim_nu" ,  lambda , nu , psi0 , psi0b , alpha0 , alpha0b , align.scale , align.pow ,
					wgt , eps , group.combis-1 , PACKAGE="sirt")	
						}