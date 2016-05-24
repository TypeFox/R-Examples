

################################################
# item fit statistics based on WLEs
msq.itemfitWLE <- function( tamobj , fitindices = NULL , ...  ){ 
	s1 <- Sys.time()
	parmlabel <- fitindices	
	res <- IRT.residuals( tamobj , ... )
	resp <- tamobj$resp
	I <- ncol(resp)
	
	#********************************
	# item indices	
	
	stand_residuals <- res$stand_residuals
	X_var <- res$X_var
	#********************************
    # calculate fit statistics
	Z <- stand_residuals^2
	dfr_data <- data.frame( "item" = colnames(resp) , 
					"N" = colSums( 1 - is.na(resp) ) )
	dfr_data$sum_var <- colSums( X_var , na.rm=TRUE)
	# outfit
	dfr_data$Outfit <- colSums( Z , na.rm=TRUE ) / dfr_data$N 
	# infit
	dfr_data$Infit <- colSums( X_var * Z  , na.rm=TRUE ) / dfr_data$sum_var
	
	#******************
	# fit statistic for subgroups
	dfr_parm <- NULL
	if ( ! is.null( parmlabel) ){
			parms <- sort( unique(parmlabel) )
			ind <- match( parmlabel , parms )
			dfr_parm <- data.frame( "parm" = parms )
			dfr_parm$N <- stats::aggregate( dfr_data$N , list(ind) , sum )[,2]
			dfr_parm$sum_var <- stats::aggregate( dfr_data$sum_var , list(ind) , sum )[,2]
			
			out1 <- stats::aggregate( dfr_data$N * dfr_data$Outfit , list(ind) , sum )[,2]
			dfr_parm$Outfit <- out1 / dfr_parm$N
			out1 <- stats::aggregate( dfr_data$sum_var * dfr_data$Infit , list(ind) , sum )[,2]
			dfr_parm$Infit <- out1 / dfr_parm$sum_var			
			dfr_parm$sum_var <- NULL
						}
	dfr_data$sum_var <- NULL	
		# summary statistics
		vars <- c("Outfit" , "Infit")
		dfr2a <- data.frame( "fit" = vars , "M" = colMeans(dfr_data[,vars]) ,
					"SD" = apply( dfr_data[,vars] , 2 , stats::sd ) )
		dfr2b <- NULL

		if ( ! is.null(dfr_parm) ){			
			vars <- c("Outfit" , "Infit")
			dfr2b <- data.frame( "fit" = vars , "M" = colMeans(dfr_parm[,vars]) ,
						"SD" = apply( dfr_parm[,vars] , 2 , stats::sd ) )
								}
		s2 <- Sys.time()
		v1 <- c(s1 , s2 )					
	res1 <- list( fit_data = dfr_data , fit_parm = dfr_parm , 
				fitindices = fitindices ,
				fit_data_summary = dfr2a , fit_parm_summary = dfr2b ,
				time = v1 )
	class(res1) <- "msq.itemfitWLE"
	return(res1)
		}
####################################################		