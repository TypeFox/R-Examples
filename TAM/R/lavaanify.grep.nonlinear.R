


#*** grep for nonlinear terms for factors and regressions
lavaanify.grep.nonlinear <- function( lavmodel , items ){	
		lavmodel <- gsub( ";" , "\n" , paste(lavmodel) )
		syn <- unlist( strsplit( paste(lavmodel) , split="\n") )
		syn <- gsub(" " , "" , syn )
		syn <- syn[ syn != "" ]
		ind <- substring( syn , 1 , 1 ) != "#"
		syn <- syn[ind]
		
		dfr <- data.frame("index" = seq(1,length(syn) ) , "syntax" = syn )		
		dfr$nonlin_loadings <- 0
		dfr$nonlin_regr <- 0
		dfr$nonlin_terms <- ""
		dfr$syntax_mod <- dfr$syntax
				
		#**********************************
		#**** grep for nonlinear factor loadings
		ind1 <- grep( "=~" , syn  )
		ind2 <- grep("I(" , syn , fixed=TRUE)
		ind <- intersect(ind1 , ind2 )
		if ( length(ind) > 0 ){
			dfr$nonlin_loadings[ind] <- 1
			l1 <- unlist(lapply( strsplit( syn[ind]  , split="=~" ) , FUN = function(hh){ hh[1] } ))
			dfr$nonlin_terms[ ind ] <- l1
						}
		#************************************				
		#***** grep for regressions on factors				
		ind1 <- grep( "~" , syn  )
		ind1a <- grep( "~~" , syn  )
		ind1b <- grep( "=~" , syn  )
		ind1 <- setdiff( ind1 , c(ind1a,ind1b) )					
		ind2 <- grep("I(" , syn , fixed=TRUE)
		ind <- intersect(ind1 , ind2 )
		if ( length(ind) > 0 ){
			dfr$nonlin_regr[ind] <- 1						
			l1 <- unlist(lapply( strsplit( syn[ind]  , split="~" ) , FUN = function(hh){ hh[2] } ))
#			l2 <- sapply( l1 , FUN = function(ii){
		    for (kk in seq(1,length(ind) ) ){
					ii <- l1[kk]	
					ii <- unlist( strsplit( ii , split="+" , fixed=TRUE)	)
					ii <- grep( "I(" , ii , fixed=TRUE , value=TRUE)
					ii2 <- strsplit( ii , split="I(" , fixed=TRUE )
					ii2 <- unlist( lapply( ii2 , FUN = function(hh){
										hh[ length(hh) ] } ) )
					ii2 <- paste0( "I(" , ii2 )
					ii2 <- paste0( ii2 , collapse= " ## " )
					dfr$nonlin_terms[ ind[kk] ] <- ii2
								}
						}
		#**** collect all nonlinear terms
		nt1 <- unlist( strsplit( paste( dfr$nonlin_terms ) , split=" ## " ) )
		nt1 <- unique( nt1 )
		N1 <- length(nt1)
		dfr2 <- NULL
		if ( N1 > 0 ){
			dfr2 <- data.frame( "nonlin_fac" = nt1 , "nonlin_fac_rename" = paste0( "NLFAC_" , 1:N1 ) )
			for ( ii in 1:N1){
			   h2 <- paste(dfr$syntax_mod)
					dfr$syntax_mod <- gsub( paste(dfr2$nonlin_fac)[ii] , paste( dfr2$nonlin_fac_rename )[ii] , 
								h2 ,  fixed=TRUE )
							}
				dfr2$exo <- 0							
					}
					
		syn <- paste0( dfr$syntax_mod , collapse = "\n" )				
		# exogeneous nonlinear terms
		if (N1 > 0 ){
			dfr2a <- dfr2[ ! ( paste(dfr2$nonlin_fac) %in% paste( dfr$nonlin_terms)[ dfr$nonlin_loadings == 1] ) , ]			
			if ( nrow(dfr2a) > 0 ){
				items <- c( items , paste(dfr2a$nonlin_fac_rename ) )				
				dfr2$exo[ paste(dfr2$nonlin_fac_rename ) %in% paste(dfr2a$nonlin_fac_rename ) ] <- 1
				}
					
					}
					
		if ( is.null(dfr2) ){ dfr <- NULL }		
		res <- list("lavmodel" = syn , "nonlin_syntable" = dfr , 
				"nonlin_factors"=dfr2 , "items" = items )		
		return(res)	
				}
