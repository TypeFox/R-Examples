
############################################
# permutation test for LSEM model
lsem.permutationTest <- function( lsem.object, B=1000 , residualize = TRUE ,
			verbose = TRUE ){

		s1 <- Sys.time()	
		CALL <- match.call()		
		
		lavaan.args <- lsem.object$lavaan.args
		entr <- c( "lavmodel" , "data" , "h" , "moderator.grid" , "moderator" , "eps" ,
              "fit_measures")
		object <- lsem.object
		arglist <- list()
		EE <- length(entr)
		for (ee in 1:EE){
			arglist[[ entr[ee] ]] <- object[[ entr[ee] ]]
						}
		arglist2 <- lsem.object$lavaan.args
		NL <- length(arglist2)
		if (NL > 0){
		for (ll in 1:NL){
			arglist[[ names(arglist2)[ll] ]] <- arglist2[[ names(arglist2)[ll] ]]
						}
					}
		arglist$residualize <- residualize
		arglist$verbose <- FALSE
		arglist$standardized <- lsem.object$standardized
		arglist$standardized_type <- lsem.object$standardized_type
	
		data1 <- data0 <- object$data		
		parameters <- object$parameters
		parameters_summary <- object$parameters_summary
		moderator <- object$moderator
		
		#******************************************
		# start permutation test
		
		parameters_permutation <- matrix( NA , nrow(parameters) , ncol=B)
		parameters_summary_M <- matrix( NA , nrow(parameters_summary) , ncol=B)
		rownames(parameters_summary_M) <- parameters_summary$par
		parameters_summary_MAD <- parameters_summary_SD <- parameters_summary_M
		parameters_summary_lin_slo <- parameters_summary_M
		
				
		if ( verbose ){ 
			cat("Permutation test LSEM \n")
						}
		
		for (bb in 1:B){
			# bb <- 1
			if (verbose){ 
				cat( bb , " ")
				if ( bb %% 20 == 0 ){ cat("\n") }
				flush.console();
						}
			data1[ , moderator ] <- sample( data0[ , moderator ] )
            arglist$data <- data1			
			res0 <- base::do.call( lsem.estimate , arglist )
			parameters_permutation[, bb] <- res0$parameters$est
			parameters_summary_M[,bb] <- res0$parameters_summary$M
			parameters_summary_SD[,bb] <- res0$parameters_summary$SD
			parameters_summary_MAD[,bb] <- res0$parameters_summary$MAD
			parameters_summary_lin_slo[,bb] <- res0$parameters_summary$lin_slo
							}
		if (verbose){ cat("\n") }					
					
		#*****************
		# create global test statistics
		teststat <- data.frame( "par" = parameters_summary$par )
		teststat$M <- parameters_summary$M
		teststat$SD <- parameters_summary$SD
		teststat$SD_p <- rowMeans( parameters_summary_SD >= parameters_summary$SD )
		teststat$MAD <- parameters_summary$MAD
		teststat$MAD_p <- rowMeans( parameters_summary_MAD >= parameters_summary$MAD )
		teststat$lin_slo <- parameters_summary$lin_slo
		p1 <- rowMeans( parameters_summary_lin_slo >= parameters_summary$lin_slo )
		p2 <- rowMeans( parameters_summary_lin_slo <= parameters_summary$lin_slo )		
		teststat$lin_slo_p <- 2*ifelse( p1 < p2 , p1 , p2 )
		teststat$lin_slo_p <- ifelse( teststat$lin_slo_p > 1 , 
									       1 , teststat$lin_slo_p )				
		#********************
		# pointwise statistics
		rownames(parameters_permutation) <- rownames(parameters)

		parameters_pointwise_test <- parameters[ , c("grid_index" ,
				                "moderator","par","parindex" ) ]
		parameters_pointwise_test$est <- parameters$est - parameters_summary$M[ parameters$parindex ]	

		par_pointwise_perm <- parameters_permutation - parameters_summary_M[ parameters$parindex , ] 

		p1 <- rowMeans( parameters_pointwise_test$est >= par_pointwise_perm )
		p2 <- rowMeans( parameters_pointwise_test$est <= par_pointwise_perm )		
		parameters_pointwise_test$p <- 2*ifelse( p1 < p2 , p1 , p2 )
		parameters_pointwise_test$p <- ifelse( parameters_pointwise_test$p > 1 , 
									       1 , parameters_pointwise_test$p )
		
		
		s2 <- Sys.time()
		
		res <- list( "teststat" = teststat , 
					"parameters_pointwise_test" = parameters_pointwise_test ,
					 "parameters" = parameters ,
					 "parameters_permutation" = parameters_permutation ,
					 "parameters_summary" = parameters_summary , 
					 "parameters_summary_M" = parameters_summary_M ,
					 "parameters_summary_SD" = parameters_summary_SD ,
					 "parameters_summary_MAD" = parameters_summary_MAD ,
					 "parameters_summary_lin_slo" = parameters_summary_lin_slo	,
					 "par_pointwise_perm" = par_pointwise_perm , 
					 "moderator.density" = object$moderator.density	,
					 "moderator" = object$moderator	,
					 "moderator.grid" = object$moderator.grid ,
					 "h" = object$h , "bw" = object$bw , "N" = object$N , 
					 "B"= B , "s1"=s1 , "s2"=s2 ,
					 "lavmodel" = object$lavmodel , CALL = CALL
							)
		class(res) <- "lsem.permutationTest"							
		return(res)
					}
################################################################					