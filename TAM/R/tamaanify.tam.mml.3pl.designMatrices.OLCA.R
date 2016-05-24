

######################################
# OLCA
tamaanify.tam.mml.3pl.designMatrices.OLCA <- function( res ){
	anlist <- res$ANALYSIS.list
	items <- colnames(res$resp)
	I <- length(items)
	itemtable <- res$items

	
		# stopping criterion for polytomous data
		maxK <- res$maxcat
		if ( maxK > 1 ){
		   stop(paste0( "Ordered latent class analysis (OLCA) is only\n" , 
		         "  available for dichotomous data!") )
						}
		TP <- ncl <- anlist$NCLASSES
		res$theta.k <- diag(TP)		
	   # design matrix E 
	   # E ( item , category , skill class , parameter )	
		Nparm <- TP*sum( itemtable$ncat - 1 ) + 1		
		E <- array( 0 , dim=c(I,maxK+1,TP , Nparm) )
		dimnames(E)[[1]] <- items
		dimnames(E)[[2]] <- paste0( "Cat" , 1:(maxK+1) )
		dimnames(E)[[3]] <- paste0("Class",1:TP)
		# E[ , -1 , , Nparm ] <- -99
				
		# define E item parameter table
		h1 <- sapply( itemtable$ncat - 1 , FUN = function(vv){
					1:vv } , simplify=FALSE)
		dfr <- data.frame( "item" = rep( itemtable$item , itemtable$ncat-1) , 
					"cat" = unlist(h1)	)
		dfr$parm <- paste0( dfr$item , "_Cat" , dfr$cat)
		dfr$indexparm <- seq( 1 , nrow(dfr) )				
		itempars <- NULL
		for (tt in 1:TP){
			# tt <- 1
			itempars <- c(itempars , paste0( dfr$parm , "_Cl" , tt ) )
						}
		itempars <- c(itempars , "one")	
		dimnames(E)[[4]] <- itempars
		M1 <- nrow(dfr)

		for (mm in 1:M1){
		    for (tt in 1:TP){
				# tt <- 1
				E[ paste(dfr$item[mm]) , dfr$cat[mm]+1  , 
							tt , mm + (1:tt-1)*M1 ] <- 1
							}
						}
		# set prior distributions
		if ( res$maxcat == 1 ){
			M1 <- nrow(dfr)
			gammaslope.prior <- matrix( c( 0 , 4 , -Inf , Inf ) , nrow=Nparm , ncol=4,
						byrow=TRUE )
			pospars <- seq(M1+1,TP*M1)
			gammaslope.prior[ pospars , 1] <- 0.5	
			gammaslope.prior[ pospars , 3] <- 0	
				}
		res$gammaslope.prior <- gammaslope.prior
																						
		# search for items which do not have the maximum number of categories
		items1 <- which( itemtable$ncat - 1 < maxK )
        for ( ii in items1 ){
			for ( hh in seq( itemtable$ncat[ii]  , maxK ) ){
				E[ ii , hh + 1, 1:TP , Nparm ] <- -99				
							}
						}

		res$E <- E
		res$notA <- TRUE		
		
		# fixed gammaslope parameter
		gammaslope.fixed <- cbind( Nparm , 1 )
		res$gammaslope.fixed <- gammaslope.fixed		
 
	 return(res)
		}
