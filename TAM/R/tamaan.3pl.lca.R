
#######################################################################
# tamaan 3PL LCA module
tamaan.3pl.lca <- function( res0 , anal.list , con , ... ){
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
		Edes <- NULL
		
		if (NSTARTS[1] > 0 ){
		for (nn in 1:(NSTARTS[1]) ){		
			gammaslope <- c( stats::qlogis( stats::runif( dim(res0$E)[4] - 1 ) ) , 1 )
			TP <- dim(res0$E)[3]
			
			# delta.inits
			if (nn==1){ delta.inits <- NULL }
			if (nn>1){
			    G <- ncol(res$delta)
				delta.inits <- inits.delta.lca( G , TP , nn  )
					}				
				
					
			res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , gammaslope=gammaslope ,  
						notA= res0$notA , control=con0 , delta.inits=delta.inits , 
						gammaslope.fixed = res0$gammaslope.fixed , 
						Edes = Edes , 	... )
			Edes <- res$Edes
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
			gammaslope <- c( stats::qlogis( stats::runif( dim(res0$E)[4] - 1 ) ) , 1 )				
			delta.inits <- NULL
							}
			
		res <- tam.mml.3pl(resp= res0$resp , E=res0$E , skillspace="discrete" , 
						theta.k= res0$theta.k , gammaslope=gammaslope ,  
						gammaslope.fixed = res0$gammaslope.fixed ,
						notA= res0$notA , delta.inits = delta.inits ,  control=con , 
						Edes = Edes , ... )
		# LCA probabilities
        res$lcaprobs <- .extract.lcaprobs(res)
		
		res$tamaan.method <- "tam.mml.3pl"
        return(res)
			}
#######################################################################

##########################################################
# extract LCA probabilities
.extract.lcaprobs <- function(res){	
		rprobs <- res$rprobs
		resp <- res$resp
		rpdim <- dim(rprobs)
		ncat <- rpdim[2]
		I <- rpdim[1]
		TP <- rpdim[3]
		obji <- NULL
		for (hh in 1:ncat){
			obji <- rbind( obji , rprobs[,hh,] )
						}
		colnames(obji) <- paste0("Class", 1:TP )
		obji <- data.frame( "item" = rep( colnames(resp) , ncat) ,
						"itemno" = rep( 1:I , ncat) , "Cat" =
						rep(1:ncat , each=I)-1 , obji )
		obji <- obji[ order( obji$itemno ) , ]		
		rownames(obji) <- NULL
		return(obji)
				}
#######################################################

######################################################
# inits delta parameters
inits.delta.lca <- function( G , TP , nn  ){
		delta.inits <- matrix( 1/TP , nrow=TP , ncol=G )
		delta.inits <- delta.inits + nn/8 * stats::runif(TP*G)
		delta.inits <- delta.inits / colSums(delta.inits)
		delta.inits <- log(delta.inits)
		return(delta.inits)
				}
########################################################				