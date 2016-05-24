
############################################
# model prior parsing
tamaanify.modelprior <- function(res){
	t1 <- res$tammodel.dfr
	gammaslope.prior <- NULL	
	xsi.prior <- NULL
	guess.prior <- NULL
	ind <- which( t1$syn == "MODELPRIOR:" )
	dfr <- matrix( 0 , nrow=0, ncol=2 )

	#********************

	if ( length(ind) > 0 ){
		n1 <- t1[ ind , "part_begin" ]
		t1 <- t1[ ( t1$part_begin == n1 ) , ]
		t1 <- t1[-1,]
		syn <- paste(t1$syn)
		l1 <- strsplit( syn , split="~" , fixed=TRUE )
		dfr <- data.frame( "parm" = unlist(lapply( l1 , FUN = function(uu){ uu[1] }  ) ) )
		dfr$prior <- unlist(lapply( l1 , 
				FUN = function(uu){ uu[2] }  ) ) 
		prior <- paste(dfr$prior)
		prior <- gsub("(" , "#" , prior , fixed=TRUE )
		prior <- gsub(")" , "#" , prior , fixed=TRUE )
		prior <- gsub("," , "#" , prior , fixed=TRUE )
		l1 <- strsplit( prior , split="#" , fixed=TRUE )
		dfr$dist <- ""
		dfr$par1 <- NA
		dfr$par2 <- NA
		dfr$par3 <- NA
		dfr$par4 <- NA
		dfr$dist <- unlist( lapply( l1 , FUN = function(uu){ uu[1] } ) )
		dfr$par1 <- as.numeric(unlist( lapply( l1 , FUN = function(uu){ uu[2] } ) ))		
		dfr$par2 <- as.numeric(unlist( lapply( l1 , FUN = function(uu){ uu[3] } ) ))		
		dfr$par3 <- as.numeric(unlist( lapply( l1 , FUN = function(uu){ 
					if ( length(uu) >= 4 ){ gg <- uu[4] } else {gg <- NA}
					return(gg) } ) ))		
		dfr$par4 <- as.numeric(unlist( lapply( l1 , FUN = function(uu){ 
					if ( length(uu) >= 5 ){ gg <- uu[5] } else {gg <- NA}
					return(gg) } ) ))


					
		#********
		# type of prior parameter
		E <- res$E
		Epars <- dimnames(E)[[4]]		
		dfr$gammaslope_index <- match( paste(dfr$parm) , Epars )
		NP <- length(Epars)
		gammaslope.prior <- matrix( 0 , nrow=NP , ncol=2 )
		gammaslope.prior[,2] <- 100
		rownames(gammaslope.prior) <- Epars
		dfr0 <- dfr[ ! is.na( dfr$gammaslope_index ) , ]
		if ( nrow(dfr0) > 0 ){
			gammaslope.prior[ dfr0$gammaslope_index , ] <- 
					as.matrix( dfr0[ , c("par1","par2") ] )
							}	
		if ( nrow(gammaslope.prior) == 0 ){
				gammaslope.prior <- NULL
									}							
		#*** index A parameters
		A <- res$A
		dfr$A_index <- match( paste(dfr$parm) , dimnames(A)[[3]] )
		Aparm <- dimnames(A)[[3]]
		NAparm <- length(Aparm)
		xsi.prior <- matrix( NA , nrow=NAparm , ncol= 2 )
		xsi.prior[,1] <- 0
		xsi.prior[,2] <- 1000
		rownames(xsi.prior) <- Aparm
		dfr0 <- dfr[ ! is.na( dfr$A_index ) , ]
		if ( nrow(dfr0) > 0 ){
			xsi.prior[ dfr0$A_index , ] <- 
					as.matrix( dfr0[ , c("par1","par2") ] )		
							}		
					}
	
		#**** index priors guessing parameters			
		if ( nrow(dfr) > 0 ){
			guess.prior <- NULL
			lavpartable <- res$lavpartable
			lav1 <- lavpartable[ lavpartable$op == "?=" , ]
			dfr$guess_index <- match( paste(dfr$parm) , paste( lav1$label ) )
						 
	
		if ( sum( 1 - is.na( dfr$guess_index )  ) > 0 ){
			items <- colnames(res$resp)
			I <- length(items)
			guess.prior <- matrix( 1 , nrow=I , ncol=2 )
			rownames(guess.prior) <- items
			est.guess <- res$est.guess
			for (kk in 1:2){
				guess.prior[,kk] <- ifelse( est.guess == 0 , 0 , guess.prior[,kk] )
							}
			
			dfr0 <- dfr[ ! is.na( dfr$guess_index ) , ]
			NP <- nrow(lav1)
			lav1 <- lav1[ paste(lav1$label) %in% paste(dfr0$parm) , ]
			ind <- match( paste(lav1$label) , paste(dfr0$parm)   )
			dfr11 <- dfr0[ ind , ]
			rownames(dfr11) <- paste0( lav1$lhs )
			guess.prior[ rownames(dfr11) , ] <- as.matrix(dfr11[ , c("par1" , "par2") ])
								}
							}
			

		if ( is.null( res$gammaslope.prior) ){
			res$gammaslope.prior <- gammaslope.prior
											}
		res$xsi.prior <- xsi.prior
		res$guess.prior <- guess.prior
		return(res)
		}