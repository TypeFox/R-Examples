
#######################################################################
# tamaan 3PL mixture module
tamaan.3pl.mixture <- function( res0 , anal.list , con , ... ){
		if ( ! is.null( anal.list$NSTARTS ) ){
				NSTARTS <- anal.list$NSTARTS 
						} else { 
				NSTARTS <- c(0,0) 
						}
		
		
		#*** initial gammaslope estimate
		
		# different starts if NSTARTS > 0
		con0 <- con
		con0$maxiter <- NSTARTS[2]
		con0$progress <- FALSE
		devmin <- 1E100
		
		itempartable <- res0$itempartable_MIXTURE
		itempartable.int <- itempartable[ itempartable$int == 1 , ]
		itempartable.slo <- itempartable[ itempartable$slo == 1 , ]
		
		gammaslope0 <- itempartable$val
		
		resp <- res0$resp
		items0 <- res0$items
		
		# initial values
		I <- ncol(resp)
		beta0 <- sapply( 1:I , FUN = function(ii){
			# ii <- 2
			ncat.ii <- items0[ii , "ncat"] - 1
			l1 <- rep(0,ncat.ii)
			for (hh in 1:ncat.ii){
				l1[hh] <- stats::qlogis( mean( resp[,ii] >= hh , na.rm=TRUE ) / ncat.ii )
								}
			return(l1)
					} )
		beta0 <- unlist( beta0)		
		B0 <- length(beta0)		
	    ncl <- anal.list$NCLASSES
		
		if (NSTARTS[1] > 0 ){
		   for (nn in 1:(NSTARTS[1]) ){
		
			gammaslope <- gammaslope0			
			gammaslope[ itempartable.int$index ] <-	rep( beta0 , ncl ) + 
							stats::rnorm( ncl*B0 , mean=0, sd = log(1+nn^(1/5) ) )
			
			N0 <- nrow(itempartable.slo)
			if ( ! res0$raschtype ){
				gammaslope[ itempartable.slo$index ] <-	stats::runif( N0 , max(.2,1-nn/5) , min( 1.8 , 1+nn/5) )			
								}
			
			
			# delta.inits
			if (nn==1){ delta.inits <- NULL }
#			if (nn>1){
#			    G <- ncol(res$delta)
#				delta.inits <- inits.delta.lca( G , TP , nn  )
#					}				


			res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , gammaslope=gammaslope ,  
						gammaslope.constr.V = res0$gammaslope.constr.V,
						gammaslope.constr.c = res0$gammaslope.constr.c,						
						notA= TRUE , control=con0 , delta.inits=delta.inits , 
						delta.designmatrix = res0$delta.designmatrix ,
						delta.fixed = res0$delta.fixed , 
						gammaslope.fixed = res0$gammaslope.fixed , 
						... )
            if (con$progress){
				cat( paste0( "*** Random Start " , nn  ,
						" | Deviance = " , round( res$deviance , 2 ) , "\n") )
				utils::flush.console()
						}
			if ( res$deviance < devmin ){
				devmin <- res$deviance
				gammaslope.min <- res$gammaslope
				delta.min <- res$delta
						}
					}
				}
		
		# use inits or best solution from random starts		
		if (NSTARTS[1] > 0 ){
			gammaslope <- gammaslope.min
			delta.inits <- delta.min
					} else {
			gammaslope <- NULL				
			delta.inits <- NULL
							}
		
		res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , gammaslope=gammaslope ,  
						gammaslope.fixed = res0$gammaslope.fixed ,
						gammaslope.constr.V = res0$gammaslope.constr.V,
						gammaslope.constr.c = res0$gammaslope.constr.c,												
						notA= TRUE , delta.inits = delta.inits ,  
						delta.fixed = res0$delta.fixed , 
						control=con , 
						delta.designmatrix = res0$delta.designmatrix ,
						... )
		
		#*****************************************
		# processing output
				
		# probabilities mixture distributions
        itempartable <- res0$itempartable_MIXTURE
		theta_MIXTURE <- res0$theta_MIXTURE
		TG <- nrow(theta_MIXTURE)
		TP <- ncl*TG
		pi.k <- res$pi.k
		D <- ncol(theta_MIXTURE )
		G <- 1
		
		# mixture probabilities
		probs_MIXTURE <- rep(NA,ncl)
		names(probs_MIXTURE) <- paste0("Cl" , 1:ncl )
		moments_MIXTURE <- as.list( 1:ncl )
		for (cl in 1:ncl){
			cl.index <- 1:TG + (cl-1)*TG
			probs_MIXTURE[cl] <- sum(pi.k[ cl.index , 1 ] )
			pi.ktemp <- pi.k[ cl.index ,,drop=FALSE]
			pi.ktemp <- pi.ktemp / colSums( pi.ktemp)
			moments_MIXTURE[[cl]] <- .mml.3pl.distributionmoments( D =D , 
					G =G , pi.k=pi.ktemp , theta.k=theta_MIXTURE )
					}
		# item parameters
		res$probs_MIXTURE <- probs_MIXTURE
		res$moments_MIXTURE <- moments_MIXTURE
		
		ipar <- res0$itempartable_MIXTURE
		p11 <- strsplit( paste(ipar$parm) , split="_Cl" )
		ipar$parm0 <- unlist(  lapply( p11 , FUN = function(pp){ pp[1] } ) )
		
		ipar$est <- gammaslope[ ipar$index ]
		# res$itempartable1_MIXTURE <- ipar
		res$gammaslope <- gammaslope
		
		# second item parameter table
		ipar2 <- ipar[ ipar$Class == 1 , c("item" , "parm0")]
		colnames(ipar2)[2] <- "parm"
		for (cl in 1:ncl){
			#cl <- 1		
			ipar2[ , paste0("Cl" , cl ) ] <- ipar[ ipar$Class == cl , "est" ]
						}
		res$itempartable_MIXTURE <- ipar2

		res$tamaan.method <- "tam.mml.3pl"
        return(res)
			}
#######################################################################
