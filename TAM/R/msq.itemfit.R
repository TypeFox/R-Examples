
#######################################
# Item fit mean squares statistics
msq.itemfit <- function( object , fitindices=NULL ){
		s1 <- Sys.time()
		#--- collect necessary input
		resp <- IRT.data(object)
# a0 <- Sys.time()		
		res <- predict(object)
		irf1 <- IRT.irfprob(object)
		irf1[ is.na(irf1) ] <- 0
# cat("predict function\n"); a1 <- Sys.time() ; print(a1-a0); a0 <- a1		
		post1 <- IRT.posterior(object)
		NI <- dim(res$probs.categ)
		I <- NI[4]
		TP <- NI[3]
		K <- NI[2]
		N <- NI[1]
		fititem <- FALSE
		if ( is.null(fitindices) ){
			fitindices <- seq(1,I)
			fititem <- TRUE
								}		
#		fitgroups <- unique( fitindices)
#		fitindices <- match( fitindices , fitgroups )
#		FF <- length( fitgroups )		

					
		#---- fit statistics in R		
# cat("Rest\n"); a1 <- Sys.time() ; print(a1-a0); a0 <- a1		
#		dfr <- msq.itemfit.R( dfr , FF , fitindices , fitgroups ,
#					res , post1 , N , TP , I , K , resp )		
# cat("msq in R\n"); a1 <- Sys.time() ; print(a1-a0); a0 <- a1				
	
		#---- include fit statistic in Rcpp
		
		fitgroups <- unique( fitindices)
		FF <- length(fitgroups)
		fitindices <- match( fitindices , fitgroups )
		fitIndexM0 <- fitIndexM <- cbind( seq(1 , length(fitindices) ) , fitindices )
		m1 <- order(fitIndexM[,2] )
		fitIndexM <- as.matrix( fitIndexM[ m1 , ] - 1) 
		fitIndexM0 <- as.matrix( fitIndexM0[ m1 , ] ) 
		fitIndexTable <- matrix( 0 , nrow=FF , ncol=3 )
		fitIndexTable[,1] <- 1:FF
		fitIndexTable[,2] <- sapply( 1:FF , FUN = function(ff){ 
					g1 <- which( fitIndexM0[,2] == ff ) 
					return(g1[1])
							} )
		fitIndexTable[,3] <- sapply( 1:FF , FUN = function(ff){ 
					g1 <- max( which( fitIndexM0[,2] == ff )  )
					return(g1)
							} ) 
		FIT_ <- as.matrix( fitIndexTable - 1 )

		dfr <- data.frame( "fitgroup" = fitgroups )		
		if (fititem){
			dfr <- cbind( "item"= colnames(resp) , dfr )
					}
		
		
		irf1_ <- as.numeric(irf1)
		res0 <- .Call( "msq_itemfit" , resp , irf1_ , K , TP , post1 , FIT_ , fitIndexM ,
					               PACKAGE="TAM")$dfr_fit
		res0 <- as.data.frame(res0)
		colnames(res0) <- c("Outfit" , "Outfit_t" , "Infit" , "Infit_t")
		dfr <- cbind( dfr , res0)
		
		#---   compute p values
		dfr$Outfit_p <- 2 * stats::pnorm( -abs( dfr$Outfit_t )	)			
		dfr$Infit_p <- 2 * stats::pnorm( -abs( dfr$Infit_t ))	
		# arrange columns in data frame
		cdfr <- colnames(dfr)
		ind <- c( grep( "Outfit" , cdfr ) ,  grep( "Infit" , cdfr ) )
		ind <- c( setdiff( seq(1,ncol(dfr) ) , ind ) , ind )
		dfr <- dfr[,ind]	
		# summary statistic
		vars <- c("Outfit" , "Infit")
		dfr2 <- data.frame( "fit" = vars , "M" = colMeans(dfr[,vars]) ,
					"SD" = apply( dfr[,vars] , 2 , stats::sd ) )
		s2 <- Sys.time()
		v1 <- c(s1 , s2 )
		res <- list( "itemfit" = dfr , "summary_itemfit"=dfr2 ,
					time=v1 )
		class(res) <- "msq.itemfit"
		return(res)			
					}
######################################################					

##########################################################
# fit statistic in R
msq.itemfit.R <- function( dfr , FF , fitindices , fitgroups ,
					res , post1 , N , TP , I , K , resp ){

		#----------------------------------
		# loop over fit groups
		for (ff in 1:FF){			
			# ff <- 1
			ind.ff <- which( fitindices == ff )

			#*********
			# Outfit statistic

			#--- compute fit statistic
			fit0 <- rep(0,N)
			for (ii in ind.ff){
				stand.resid <- res$stand.resid[,,ii]
				fit1 <- rowSums( post1 * stand.resid^2 , na.rm=TRUE)
				fit0 <- fit1 + fit0
							 }
			fit0 <- sum(fit0) / sum( 1 - is.na( resp[,ind.ff] ) )
			dfr[ff,"Outfit"] <- fit0

			#--- compute inference

			v1 <- rep(0,N)
			for (ii in ind.ff){
				probs.ii <- res$probs.categ[,,,ii]
				kurt.ii <- array( 0 , dim=c(N,TP) )
				exp.ii <- res$expected[,,ii]
				for (kk in 1:K){
				kurt.ii <- kurt.ii + probs.ii[,kk,] * ( kk - 1 - exp.ii )^4
								}
				
				v0 <- kurt.ii / res$variance[,,ii]^2 
				v0 <- rowSums( post1 * v0 , na.rm=TRUE)
				v1 <- v1 + v0
							}
			N1 <- sum( 1-is.na(resp[,ind.ff] ) )
			qi <- sum( v1 / N1^2 ) - 1/N1
	
			# this seems to be the adequate formula
			dfr[ff,"Outfit_t"] <- ( fit0^(1/3)-1 )* 3 / sqrt(qi) + sqrt(qi) / 3 

			#***********
			# Infit statistic
			term1 <- term2 <- rep(0,N)
			for (ii in ind.ff){
				stand.resid <- res$stand.resid[,,ii]
				variance <- res$variance[,,ii]
				term1 <- term1 + rowSums( post1 * stand.resid^2 * variance , na.rm=TRUE )
				term2 <- term2 + rowSums( post1 * variance , na.rm=TRUE )
							   }
			fit1 <- sum( term1 , na.rm=TRUE) / sum(term2 , na.rm=TRUE)
			dfr[ff,"Infit"] <- fit1
			
			
			v1 <- rep(0,N)
			for (ii in ind.ff){
				variance <- res$variance[,,ii]
				v0 <- ( kurt.ii - res$variance[,,ii]^2 )
				v0 <- rowSums( post1 * v0 )
				v1 <- v1+v0
								}
			qi <- sum(v1 , na.rm=TRUE) / ( sum( term2 , na.rm=TRUE) )^2
			# this seems to be a an adequate formula
			dfr[ff,"Infit_t"] <- ( fit1^(1/3)-1 )* 3 / sqrt(qi) + sqrt(qi) / 3 
					}	
		return(dfr)
			}
#######################################################################			
							