

######################################
# LOCLCA
tamaanify.tam.mml.3pl.designMatrices.LOCLCA <- function( res ){
	anlist <- res$ANALYSIS.list
	items <- colnames(res$resp)
	I <- length(items)
	itemtable <- res$items


		TP <- ncl <- anlist$NCLASSES
		res$theta.k <- diag(TP)
	
	   # design matrix E 
	   # E ( item , category , skill class , parameter )	
		Q <- res$Q	
		
		D <- ncol(Q)   # number of dimensions
		
		#***** number of parameters to be estimated
		lavpartable <- res$lavpartable
		
		Nparm <- TP*D + sum( itemtable$ncat - 1 ) + 1		
		maxK <- res$maxcat
		
		# define E item parameter table
		h1 <- sapply( itemtable$ncat - 1 , FUN = function(vv){
					1:vv } , simplify=FALSE)
		dfr <- data.frame( "item" = rep( itemtable$item , itemtable$ncat-1) , 
					"cat" = unlist(h1)	)
		dfr$parm <- paste0( dfr$item , "_Cat" , dfr$cat)
		dfr$indexparm <- seq( 1 , nrow(dfr) )				
		lav0 <- lav1 <- res$lavpartable
		ind <- which( ( paste(lav1$op  ) == "|" ) &
						( lav1$free == 0 ) )
		lav1 <- lav1[ind, ]
		dfr$dim <- unlist(apply( Q , 1 , FUN = function(vv){
					which( vv > 0 )[1] } ) )	
		dfr$fixed <- 0
		dfr$value <- NA
		
		lav0$label2 <- paste0( lav0$lhs , "_Cat" , substring(lav0$rhs , 2 ))
		lav0[ lav0$op != "|" , "label2" ] <- ""
		res$lavpartable <- lav0
		i1 <- match( paste(dfr$parm) , paste(lav0$label2) )
		dfr[ , "parm" ] <- paste(lav0[ i1 , "label" ] )
		if (nrow(lav1) > 0 ){
			dfr3 <- data.frame( "item" = paste(lav1$lhs) , 
								"cat" = substring( paste(lav1$rhs),2) ,
								"val" = lav1$ustart )
			dfr3$parm <- paste0( dfr3$item , "_Cat" , dfr3$cat )
#			dfr3$parmid <- match( dfr3$parm , dimnames(E)[[4]] )
			ind2 <- match( paste(dfr3$parm) , paste(dfr$parm) )
			dfr$fixed[ind2] <- 1
			dfr$value[ind2] <- dfr3$val 			
						}

		# check for all dimensions
		D <- max( dfr$dim )
		for (dd in 1:D){
			ind.dd <- which( dfr$dim == dd )
			if ( sum( dfr[ind.dd, "fixed"] ) == 0 ){
				dfr[ ind.dd[1] , "fixed"] <- 1
				dfr[ ind.dd[1] , "value" ] <- 0
										}
						}
		
		# itempars <- NULL
		# class locations
		dfr2 <- data.frame( "Class" = rep(1:TP ,each=D) ,
						"Dim" = rep(1:D,TP) )
		dfr2$parm <- paste0("Cl" , dfr2$Class , "_Dim" , dfr2$Dim )
		res$loclca_LOC <- dfr2
		res$loclca_ITEMS <- dfr

#		Nparm <- length(unique(dfr$parm)) + nrow(dfr2) + 1 
		Nparm <- nrow(dfr2) + 1 
				
		E <- array( 0 , dim=c(I,maxK+1,TP , Nparm) )
		dimnames(E)[[1]] <- items
		dimnames(E)[[2]] <- paste0( "Cat" , 1:(maxK+1) )
		dimnames(E)[[3]] <- paste0("Class",1:TP)		
		dimnames(E)[[4]] <- c( dfr2$parm , "one" )		
		
		M1 <- nrow(dfr)
		for (mm in 1:M1){
#		    E[ paste(dfr$item[mm]) , dfr$cat[mm]+1 , , 	paste(dfr$parm)[mm] ] <- - 1
			for (tt in 1:TP){
				E[ paste(dfr$item[mm]) , dfr$cat[mm]+1 , tt , 	
						paste0( "Cl", tt , "_Dim" ,   dfr$dim[mm] ) ] <- dfr$cat[mm]				
							}
					}

					
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
		# extact constraints
			dfr3 <- dfr[ dfr$fixed == 1 , ]
			dfr3$parmid <- match( dfr3$parm , dimnames(E)[[4]] )
			
			gammaslope.fixed1 <- as.matrix(dfr3[ , c("parmid" , "value") ])
			colnames(gammaslope.fixed1) <- colnames(gammaslope.fixed)
			gammaslope.fixed <- rbind( gammaslope.fixed , gammaslope.fixed1 )		
			gammaslope.fixed <- gammaslope.fixed[ ! is.na(gammaslope.fixed[,1] ),
									, drop=FALSE ]
			res$gammaslope.fixed <- gammaslope.fixed
			
 
	 return(res)
		}
