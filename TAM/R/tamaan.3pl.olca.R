


#######################################################################
# tamaan 3PL OLCA module
tamaan.3pl.olca <- function( res0 , anal.list , con , ... ){
		if ( ! is.null( anal.list$NSTARTS ) ){
				NSTARTS <- anal.list$NSTARTS 
						} else { 
				NSTARTS <- c(0,0) 
						}
						
		# different starts if NSTARTS > 0
		con0 <- con
		con0$maxiter <- NSTARTS[2]
		con0$progress <- FALSE
		devmin <- 1E100
		gammaslope.prior <- res0$gammaslope.prior
		Nparm <- dim(res0$E)[4]
		if (NSTARTS[1] > 0 ){
		for (nn in 1:(NSTARTS[1]) ){		
			# gammaslope <- c( qlogis( runif( dim(res0$E)[4] - 1 ) ) , 1 )				
			gammaslope <- .inits.olca.1dim(res0 , anal.list )
			
			# delta.inits
			if (nn==1){ delta.inits <- NULL }
			if (nn>1){
			    G <- ncol(res$delta)
				TP <- anal.list$NCLASSES
				delta.inits <- inits.delta.lca( G , TP , nn  )
					}						
			
			gammaslope[Nparm] <- 1
			res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , xsi.prior = res0$xsi.prior ,
						gammaslope=gammaslope ,  
						gammaslope.prior = res0$gammaslope.prior ,
						notA= res0$notA , control=con0 , 
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
	    #********************************************************
		# use inits or best solution from random starts		
		if (NSTARTS[1] > 0 ){
			gammaslope <- gammaslope.min
			delta.inits <- delta.min
				
					} else {
			gammaslope <- .inits.olca.1dim(res0 , anal.list )
			delta.inits <- NULL
							}		

				
		res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , xsi.prior = res0$xsi.prior , 
						gammaslope=gammaslope ,  
						gammaslope.prior = res0$gammaslope.prior ,
						notA= res0$notA , delta.inits = delta.inits ,  control=con , 
						... )
		# LCA probabilities
        res$lcaprobs <- .extract.lcaprobs(res)
		
		res$tamaan.method <- "tam.mml.3pl"
        return(res)
			}
#######################################################################

###########################################
# initial values 1dim OLCA
.inits.olca.1dim <- function(res0 , anal.list ){
			qm <- - stats::qlogis( colMeans( res0$resp , na.rm=TRUE) )
			TP <- anal.list$NCLASSES
			theta1 <- stats::qnorm( ( 1:TP - .5 ) / TP  )
			pr1 <- outer( theta1 , qm , "-" )
			Nparm <- dim(res0$E)[4]
			gammaslope <- rep(1,Nparm)
			N1 <- (Nparm-1)/TP
			gammaslope[ 1:N1 ] <- pr1[1,]
			for (tt in 2:TP){
			    v1 <- pr1[tt,] - pr1[tt-1,]
				L1 <- length(v1)
				gammaslope[ 1:N1+(tt-1)*N1 ] <- stats::runif(L1,0.7,1.3)*v1
					}
			return(gammaslope)
			}
######################################################
