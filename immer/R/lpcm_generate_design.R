
################################################
# design matrix LPCM
lpcm_generate_design <- function( pars_info , irtmodel , W , 
			b_const , normalization , I , maxK , nullcats ){
	
		W0 <- W
		b0 <- b_const
		n1 <- nrow(pars_info)

		if ( is.null(irtmodel) ){
			irtmodel <- "PCM"
							}
		
		if ( is.null( b_const ) ){
				b_const <- rep( 0 , n1 )
								}		
		
		
		pars_info$estpar <- 1	
		if ( nullcats == "zeroprob" ){
			pars_info$estpar <- 1*(pars_info$Freq > 0	)
			b_const[ pars_info$estpar == 0 ] <- 99
								}
								
		#*********************
		if ( nullcats != "zeroprob" ){
		    n2 <- n1
			W <- matrix( 0 , nrow=n1 , ncol=n2-1)
			rownames(W) <- paste0( pars_info$item , "_Cat" , pars_info$cat )
							}

		if ( nullcats == "zeroprob" ){
		    n2 <- sum(pars_info$estpar)
			W <- matrix( 0 , nrow=n1 , ncol=n2-1)
			rownames(W) <- paste0( pars_info$item , "_Cat" , pars_info$cat )
			if ( is.null( b_const ) ){
				b_const[ pars_info$estpar == 0 ] <- 99
								}
			irtmodel <- "PCM"
							}
							
		#------------------------
		# irtmodel == "PCM"
		if (irtmodel == "PCM"){

		   pinfo2 <- pars_info[ pars_info$estpar == 1 , ]
           n1 <- nrow(pinfo2)
    	   n2 <- n1	
		   index <- pinfo2$index
		
			# PCM: normalization = "first"
			if ( normalization == "first" ){
			   W[ cbind( index[2:n1] , 1:(n2-1) ) ] <- 1
			   colnames(W) <- rownames(W)[index[-1]	]
										}
											
			# PCM: normalization = "sum"
			if ( normalization == "sum" ){
				W[ cbind( index[1:(n1-1)] , 1:(n1-1) ) ] <- 1
				W[ index[n1] , ] <- -1 
				colnames(W) <- rownames(W)[index[-n1]]	
											}
											
								}
		
		#--------------------------
		# irtmodel == "PCM2"
		if (irtmodel == "PCM2" ){
		    items <- unique( paste(pars_info$item))
			I <- max( pars_info$itemid )
			colnames(W) <- paste0("w",1:(n1-1))
			#--- normalization == "first"
			if ( normalization == "first"){
				# items
				p1 <- pars_info[ pars_info$itemid > 1	, ]
				W[ cbind( p1$index , p1$itemid - 1 ) ] <- p1$cat
				colnames(W)[ seq(1 , I-1 ) ] <- items[-1]
											}
			#--- normalization == "sum"
			if ( normalization == "sum"){
				# items
				p1 <- pars_info[ pars_info$itemid < I	, ]
				W[ cbind( p1$index , p1$itemid  ) ] <- p1$cat
 		        colnames(W)[ seq(1 , I-1 ) ] <- items[-I]				
				p1b <- pars_info[ pars_info$itemid == I , , drop=FALSE ]
				for ( kk in seq(1,nrow(p1b) ) ){
					W[ p1b$index[kk] , 1:(I-1) ] <- - p1b$cat[kk]
												}
				
									}				
				vv <- I 
				p2 <- pars_info
				p2$param <- 0
				p2$param[ p2$cat < p2$maxK ] <- 1
				p2$param <- ( p2$param > 0 ) * ( cumsum( p2$param ) + ( vv - 1 ) )
				W[ cbind( p2$index , p2$param ) ] <- 1
				p2a <- p2[ p2$param > 0 , ]
				colnames(W)[ p2a$param ] <- paste0( p2a$item , "_Step" , p2a$cat )
							
								}

		if ( ! is.null(W0) ){
				W <- W0 
							}
		if ( ! is.null(b0) ){
				b_const <- b0
							}			
								
								
		if ( is.null( colnames(W) ) ){
			colnames(W) <- paste0("par" , 1:ncol(W) )
									}
				
		#*********************
		# output
		res <- list(W=W , b_const=b_const, irtmodel=irtmodel)
		return(res)		
					}
###############################################