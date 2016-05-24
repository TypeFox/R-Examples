
#######################################################
# collects all item parameters from mirt output in
# a data frame
mirt.wrapper.coef <- function( mirt.obj ){
    coefmirt <- coef(mirt.obj)
    items <- colnames( mirt.obj@Data$data )
    # collect item parameters from mirt output
    itempars <- NULL
	#*******************
	# single group analysis
	if ( class(mirt.obj) != "MultipleGroupClass" ){		
       coefmirt <- mirt.obj@pars   	   
	   I <- length(items)
		for (ii in 1:I){
			cii <- coefmirt[[ii]]
			cii <- cii@par   
			itempars <- unique( c(itempars , names(cii) ) )
						}
		itempars <- sort(itempars)
		IP <- length(itempars)
	
		dfr <- matrix( NA , nrow=I , ncol=IP) 
		colnames(dfr) <- itempars
		dfr <- as.data.frame(dfr)
		dfr <- data.frame("item" = items , dfr )
		rownames(dfr) <- items
		for (ii in 1:I){
			cii <- coefmirt[[ii]]
			cii <- cii@par			
			dfr[ii , names(cii) ] <- cii
						}         
		rownames(dfr) <- NULL
		if( sum(colnames(dfr) %in% "g") > 0 ){
		     dfr$g <- stats::plogis(dfr$g) 
						}
		if( sum(colnames(dfr) %in% "u") > 0 ){
		     dfr$u <- stats::plogis(dfr$u) 
						}		
				
#		coefmirt <- coef(mirt.obj)
		h1 <- mirt.obj@pars		
		h1 <- h1[[I+1]]
		h2 <- h1@par
		names(h2) <- names(h1@est)
		res <- list( "coef" = dfr , "GroupPars" = h2 )						
			}			
	#*************************
    # multiple group analysis 
    if ( class(mirt.obj) == "MultipleGroupClass" ){			
		# extract group names
		groups <- paste(mirt.obj@Data$groupNames)
		G <- length(groups)
		dfr0 <- NULL
		for (gg in groups){
		#	gg <- groups[1]
			itempars <- NULL
			for (ii in items){
				 cii <- coefmirt[[gg]][[ii]]    
				 itempars <- unique( c(itempars , colnames(cii) ) )
							}											 
			itempars <- sort(itempars)
			I <- length(items)                
			IP <- length(itempars)
			dfr <- matrix( NA , nrow=I , ncol=IP) 
			colnames(dfr) <- itempars
			dfr <- as.data.frame(dfr)
			rownames(dfr) <- items
			for (ii in items){
				   cii <- coefmirt[[gg]][[ii]]    
				   dfr[ii , colnames(cii) ] <- cii
							 }  
			colnames(dfr) <- paste0( colnames(dfr) , "." , gg )
			if (gg==groups[1] ){ dfr0 <- dfr } else { dfr0 <- cbind( dfr0 , dfr ) }
							}
		dfr <- data.frame( "item" = items , dfr0 )
		rownames(dfr) <- NULL
		coefmirt <- coef(mirt.obj)
		cm1 <- as.list(1:G)
		names(cm1) <- groups
		for (gg in groups){
			cm1[[gg]] <- coefmirt[[gg]][["GroupPars"]]
						}
		res <- list( "coef" = dfr , "GroupPars" = cm1 )		
					}		
	#*************************
    return(res)
    }
#######################################################    