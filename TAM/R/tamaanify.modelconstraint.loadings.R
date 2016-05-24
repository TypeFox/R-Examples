
##########################################################
# model constraint loadings
tamaanify.modelconstraints.loadings <- function(res){
	mdfr <- res$MODELCONSTRAINT.dfr	
	lav1 <- res$lavpartable
	Q <- res$Q
	#*** create further model constraints if there are equalities
	lav1 <- lav1[ grep( "=~" , paste(lav1$fullsyn) , fixed=TRUE) , ]
	lav1$label_ren <- paste0( lav1$rhs , "_" , lav1$lhs , "_load")
	
	labels1 <- paste(lav1$label)
	t1 <- table(labels1)
	t1 <- t1[ t1 > 1]
	labels1 <- sort( names(t1) )
	AL <- length(labels1)
	if (AL>0){
	for (aa in 1:AL){	
#		aa <- 1
		lab.aa <- labels1[aa]
		ind.aa <- which( paste(lav1$label) == lab.aa ) 
		dfr.aa <- data.frame( "index" = 999 , 
				"syn" = paste0( lav1$label_ren[ind.aa] , "==" , lab.aa ) )
		dfr.aa$derived <- 1			
		dfr.aa$trafopar <- lav1$label_ren[ind.aa]	
		dfr.aa$expanded <- 0
		dfr.aa <- cbind( dfr.aa , lav1[ ind.aa , c("lhs" , "op" , "rhs" , "fullsyn") ] )
		mdfr <- rbind( mdfr , dfr.aa )
		lav1[ ind.aa , "label" ] <- lav1[ ind.aa , "label_ren" ]		
						}
				}
	

	if ( ! is.null( mdfr ) ){
		mdfr <- mdfr[ grep( "=~" , paste(mdfr$fullsyn) , fixed=TRUE) , ]
		lav1 <- lav1[ grep( "=~" , paste(lav1$fullsyn) , fixed=TRUE) , ]
		if ( nrow(mdfr) > 0 ){
			items0 <- colnames(res$resp)
			lav1$itemno <- match( lav1$lhs , items0 )		
			lav1 <- lav1[ order(lav1$itemno) , ]		
			syn2 <- mdfr$syn
			# add parameters		
			dfr <- tamaanify.grep.linequations( syn2 )		
			lav1a <- lav1[ match( dfr$lhsparm , paste(lav1$label ) ) , ]
			dfr <- cbind( lav1a[ , c("fullsyn" , "lhs" , "op" , "rhs" ) ] , dfr )
			dfr$user <- 1			
			dfr$ustart <- NA			
			ind <- which( ! paste( lav1$label ) %in% dfr$lhsparm )
            if ( length(ind)>0 ){ 			
				lav1 <- lav1[ ind , ]
				lav1$terms <- lav1$rhsparm <- lav1$lhsparm <- lav1$label
				lav1$fac <- 1
				lav1$parm <- lav1$label		
				lav1 <- lav1[ , colnames(dfr) ]	
			    dfr <- rbind( dfr , lav1 )				
							}
	
			rownames(dfr) <- NULL			
			res$modelconstraint.loading <- dfr			
			parms <- unique( paste(dfr$parm))			
			Nparm <- length(parms)
			I <- ncol(res$resp)
			maxcat <- res$maxcat	
			D <- ncol(Q)	
			A <- array( 0 , dim=c(I, D , Nparm) )
			dimnames(A)[[1]] <- colnames(res$resp)
			dimnames(A)[[2]] <- colnames(Q)
			dimnames(A)[[3]] <- parms		
			ND <- nrow(dfr)		
			for (dd in 1:ND){
				# dd <- 1
				dfr.dd <- dfr[dd,]		
				A[ dfr.dd$rhs , dfr.dd$lhs , dfr.dd$parm ] <- dfr.dd$fac
							}
			res$L <- A	
			res$method <- "tam.mml.3pl"
			
			#*******************************************
			# definition of E design matrix
			I <- ncol( res$resp)
			maxK <- max( res$maxcat) + 1
			mload <- res$modelconstraint.loading
			Epars <- sort( unique( paste(mload$parm) ) )
			Ngam <- length(Epars)+1
			lv1 <- res$lavpartable
			itemnames <- colnames(res$resp)
			E <- array( 0 , dim=c(I,maxK,D,Ngam) )
			dimnames(E)[[1]] <- itemnames
			dimnames(E)[[3]] <- colnames(Q)
			dimnames(E)[[4]] <- c( Epars , "one")
			E[ ,,1,] <- 0
			N1 <- nrow(mload)
			items <- res$items
			for (ii in 1:N1){
			    item.ii <- paste(mload[ii,"rhs"])
				ncat <- items[ paste(items$item) == item.ii , "ncat" ]
				E[ item.ii , 2:ncat , paste(mload[ii,"lhs"]) , 
							paste(mload[ii,"parm"]) ] <- seq(1,ncat-1)*mload[ii,"fac"]
				if ( ncat < maxK ){
#					E[ item.ii , seq( ncat+1,maxK) , paste(mload[ii,"lhs"]) , "one" ] <- - 99
					E[ item.ii , seq( ncat+1,maxK) , 1:D , "one" ] <- - 99
									}							
							}
			res$E <- E	
			# gammafixed
			gammaslope.fixed <- cbind( Ngam , 1 )
			ind <- which( ( lv1$op == "=~" ) & ( lv1$free == 0 ) )
			if ( length(ind) > 0 ){
				lv1 <- lv1[ ind , ]
				gsf <- cbind( match( paste(lv1$label) , dimnames(E)[[4]] ) , lv1$ustart )
				gammaslope.fixed <- rbind( gammaslope.fixed , gsf )
								}
			res$gammaslope.fixed <- gammaslope.fixed							
							}																		
					}										
	return(res)
	}
############################################################	
