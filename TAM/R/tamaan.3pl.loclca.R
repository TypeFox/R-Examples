
#######################################################################
# tamaan 3PL LOCLCA module
tamaan.3pl.loclca <- function( res0 , anal.list , con , ... ){
		if ( ! is.null( anal.list$NSTARTS ) ){
				NSTARTS <- anal.list$NSTARTS 
						} else { 
				NSTARTS <- c(1,20) 
						}
		
		#*** initial gammaslope estimate
		
		# different starts if NSTARTS > 0
		con0 <- con
		con0$maxiter <- NSTARTS[2]
		con0$progress <- FALSE
		devmin <- 1E100

        E <- res0$E
		Nparm <- dim(E)[4]

		#*******
		# include gammaslope.prior here
		# gammaslope.prior <- cbind( rep(0,Nparm), rep(2,Nparm) )
		gammaslope.prior <- res0$gammaslope.prior
		
		if (NSTARTS[1] > 0 ){
		for (nn in 1:(NSTARTS[1]) ){		
		
			# inits loclca
			D <- ncol(res0$Q)
			TP <- nrow(res0$theta.k)

			gammaslope <- rep(0,Nparm)
			names(gammaslope) <- dimnames(E)[[4]]
			
			N1 <- nrow(res0$loclca_LOC)
			N2 <- nrow(res0$loclca_ITEMS)
			# trait locations
			K <- TP 
			locs <- stats::qnorm( seq( 1/(2*K) , 1 - 1/(2*K) , 1/K ) )
			locs <- rep( locs , each= D )
			names(locs) <- paste(res0$loclca_LOC$parm)
			if (nn==1){ locs0 <- locs }
			gammaslope[  names(locs) ] <- locs
			# inits for item parameters
			dfr <- res0$loclca_ITEMS
			G <- 1
			startprobs <- rep( 1/TP , each=TP ) + (nn-1)/5 * stats::runif( TP )
			startprobs <- startprobs / sum( startprobs )
			delta.inits <- matrix( startprobs , ncol=G )		
			# delta.inits
			if (nn>1){
				if (nn==2){ gammaslope0 <- res$gammaslope }
				gammaslope[ 1:N1 ] <- gammaslope0[1:N1] + stats::rnorm(N1, sd = nn/8 )					
				if ( nn %% 2 == 0 ){  gammaslope0[ 1:N1 ] <- stats::rnorm(N1 , sd =2 ) }
				locs <- locs0 + stats::rnorm( TP*D , sd = sqrt(nn/5) )
				gammaslope[ names(locs) ] <- locs
					}				
			res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , gammaslope=gammaslope ,  
						A= res0$A , xsi.fixed = res0$xsi.fixed , xsi.prior=res0$xsi.prior, 
						control=con0 , delta.inits=delta.inits , 
						gammaslope.fixed = res0$gammaslope.fixed , 
						gammaslope.prior = gammaslope.prior ,
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
		#**************************		
		# use inits or best solution from random starts		
		if (NSTARTS[1] > 0 ){
			gammaslope <- gammaslope.min
			delta.inits <- delta.min
					} else {
			gammaslope <- rep(0,Nparm)
			delta.inits <- NULL
							}
			
		res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , gammaslope=gammaslope ,  
						gammaslope.fixed = res0$gammaslope.fixed ,
						A= res0$A , xsi.fixed = res0$xsi.fixed , xsi.prior=res0$xsi.prior,
						delta.inits = delta.inits ,  control=con , 
						gammaslope.prior = gammaslope.prior , 
						... )
		# LCA probabilities
        res$lcaprobs <- .extract.lcaprobs(res)
		res$tamaan.method <- "tam.mml.3pl"
		# extract cluster locations
		locs <- res$gammaslope  #  [ names(locs0) ]
		locs <- locs[ names(locs0) ]
		locs <- matrix( locs , nrow=TP , ncol=D , byrow=TRUE )
		rownames(locs) <- paste0("Cl" , 1:TP )
		colnames(locs) <- colnames(res0$Q)
		locs <- as.data.frame(locs)
		locs$prob <- res$pi.k
		res$locs <- locs
        return(res)
			}
#######################################################################


