################################################################
invariance.alignment <- function( lambda , nu , wgt=NULL ,
	align.scale=c(1,1) , align.pow=c(1,1) , eps=.01 , 
	h= .001 , max.increment=.2 , increment.factor=c(1.001,1.02,1.04,1.08) ,  
	maxiter =300 , conv=.0001 , fac.oldpar = c(0.01 , .2 , .5 , .85 )  , 
	psi0.init=NULL , alpha0.init=NULL ,
	progress=TRUE ){
	#************************************************
    s1 <- Sys.time()
	max.incr0 <- max.increment
	#****
	G <- nrow(lambda)   # number of groups
	I <- ncol(lambda)   # number of items
	if ( is.null(wgt) ){  wgt <- 1+0*nu }
	
	W1 <- dim(wgt)
	wgtM <- matrix( colSums(wgt,na.rm=TRUE)  , nrow=W1[1] , ncol=W1[2] , byrow=TRUE )	
	wgtM <- wgt / wgtM
	wgt <- G * wgtM 

	
	# missing indicator matrix: 1 - no missings
	missM <- 0.5 * ( (1-is.na(lambda))+ (1- is.na(nu)) )
	wgt <- wgt * missM
	wgt[ missM == 0 ] <- 0
	
	lambda[ missM == 0 ] <- mean( lambda , na.rm=TRUE )
	nu[ missM == 0 ] <- mean( nu , na.rm=TRUE )
	
	group.combis <- t( combinat::combn( G , 2 ) )
	group.combis <- rbind( group.combis , group.combis[,c(2,1) ] )
#	group.combis <- group.combis[ order( group.combis[,1] ) , ]


    group_combis <- group.combis

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
	fopt <- rep( 1E100, 2)
	lambda <- as.matrix(lambda)
#	eps <- .001
	wgt <- as.matrix(wgt)	
	
	
	# no lambda residuals
	eps2 <- 1E-10
	lambda_constant <- FALSE
	if( stats::sd( lambda ) < 1E-10 ){
			I <- ncol(lambda)
			G <- nrow(lambda)
			lambda <- lambda + matrix( stats::runif( I*G , -eps2 , eps2 ) , nrow=G , ncol=I)
			lambda_constant <- TRUE
							}
	
	
	
	#****************************
	# define design for optimization
    gridp <- base::expand.grid( increment.factor , fac.oldpar )
	GP <- nrow(gridp)

	#**************************************************
	# OPTIMIZATION LAMBDA	
	
	if (progress){
	    cat("* OPTIMIZATION LAMBDA\n")
	    cat( paste0( "|" , paste0( rep("*" , GP ) , collapse="") , "|\n|") )
		utils::flush.console()
				}
	psi0_init <- psi0	
	for (gp in 1:GP){	
		res0 <- inv.alignment2.lambda.alg( lambda , psi0_init , nu , align.scale ,
		          align.pow , iter ,  conv , fopt_change , wgt , eps ,
		          group_combis , fopt , h , max.increment , 
				  fac.oldpar = gridp[gp,2] , maxiter ,
				  alpha0 , group.combis , G , increment.factor=gridp[gp,1])    
	
	psi0 <- res0$psi0			
	if ( res0$fopt < fopt[1] ){
		psi0.min <- psi0	
		fopt[1] <- res0$fopt
					}
		if (progress){ cat("-") ; utils::flush.console() }
					}
		if (progress){ cat("|\n") ; utils::flush.console() }	

	#**************************************************
	# OPTIMIZATION NU		
	
	# iter <- 1
	max.incr0 -> max.increment
	psi0 <- psi0.min

	
	if (progress){
	    cat("* OPTIMIZATION NU\n")
	    cat( paste0( "|" , paste0( rep("*" , GP ) , collapse="") , "|\n|") )
		utils::flush.console()
				}
    
    alpha0_init <- alpha0	
	
	for (gp in 1:GP){	

	    res0 <- inv.alignment2.nu.alg( lambda , psi0 , nu , align.scale ,
		          align.pow , iter , conv , fopt_change , wgt , eps ,
		          group_combis , fopt , h , max.increment , fac.oldpar=gridp[gp,2] , maxiter ,
				  alpha0_init , group.combis , G , increment.factor=gridp[gp,1])

		alpha0 <- res0$alpha0
	if ( res0$fopt < fopt[2] ){
		alpha0.min <- alpha0	
		fopt[2] <- res0$fopt
					}					
		if (progress){ cat("-") ; utils::flush.console() }				
					}
    if (progress){ cat("|\n") ; utils::flush.console() }			
	
	#*****************************
	# calculate item statistics and R-squared measures
	# groupwise aligned loading
	lambda.aligned <- lambda / psi0.min
	nu.aligned <- nu - alpha0.min * lambda 
	# average aligned parameter
	itempars.aligned <- data.frame("M.lambda" = colMeans(lambda.aligned) ,
			"SD.lambda" = apply( lambda.aligned , 2 , stats::sd , na.rm=TRUE ) ,
			"M.nu" = colMeans( nu.aligned ) ,
			"SD.nu" = apply( nu.aligned , 2 , stats::sd , na.rm=TRUE ) 	
				)
	rownames(itempars.aligned) <- colnames(lambda)
	lambda.resid <- lambda.aligned - 
		matrix( itempars.aligned$M.lambda , G , I , byrow=TRUE )
	nu.resid <- nu.aligned - 	
		matrix( itempars.aligned$M.nu , G , I , byrow=TRUE )
		
	# R-squared measures
	Rsquared.invariance <- c(NA,NA)
	names(Rsquared.invariance) <- c("loadings" , "intercepts" )	
#	expl <- psi0.min * matrix( itempars.aligned$M.lambda , G , I , byrow=TRUE)
#	Rsquared.invariance["loadings"] <- 
#			1 - sum( (lambda - expl)^2 , na.rm=TRUE ) / 
#			sum( ( lambda - rowMeans(lambda) )^2 , na.rm=TRUE)
#	expl <- matrix( itempars.aligned$M.nu , G , I , 
#				byrow=TRUE) + 
#			alpha0.min*matrix( itempars.aligned$M.lambda , G , I , 
#				byrow=TRUE)
#	Rsquared.invariance["intercepts"] <- 
#			1 - sum( ( nu - expl)^2 , na.rm=TRUE ) / 
#			sum( (nu - rowMeans(nu) )^2 , na.rm=TRUE)
	expl <- psi0.min * matrix( itempars.aligned$M.lambda , G , I , byrow=TRUE)
	Rsquared.invariance["loadings"] <- 
			1 - sum( (lambda - expl)^2 , na.rm=TRUE ) / 
			sum( ( lambda  )^2 , na.rm=TRUE)
	expl <- matrix( itempars.aligned$M.nu , G , I , 
				byrow=TRUE) + 
			alpha0.min*matrix( itempars.aligned$M.lambda , G , I , 
				byrow=TRUE)
	Rsquared.invariance["intercepts"] <- 
			1 - sum( ( nu - expl)^2 , na.rm=TRUE ) / 
			sum( (nu  )^2 , na.rm=TRUE)

    # correlations aligned parameters			
	rbar <- c( ai.calc.corr(t(lambda.aligned)) , ai.calc.corr(t(nu.aligned)) )			
	es.invariance <- rbind( Rsquared.invariance ,
			sqrt( 1- Rsquared.invariance ) , rbar)
	rownames(es.invariance) <- c("R2" , "sqrtU2" , "rbar")
	if (lambda_constant){
		es.invariance[,"loadings"] <- c(1,0,1)
						}
	
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
			"fopt"= fopt ,
#			"fopt.history" = fopt.history[1:(max(Niter)-1) , ] ,
			"align.scale"=align.scale , "align.pow"=align.pow ,
			"s1"=s1 , "s2"=s2)
	class(res) <- "invariance.alignment"
	return(res)
		}
#####################################
