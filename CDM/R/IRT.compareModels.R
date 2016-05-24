
###############################################################################
# compare models based on likelihood and information criteria
IRT.compareModels <- function( object , ... ){
		cl <- match.call()
		cl1 <- paste(cl)[-c(1)]
		object_list <- list(...)		
		
		#****************************
		# information criteria
		irtmodfit <- FALSE
		cm1 <- NULL
		if ( length( grep("IRT.modelfit" , class(object) ) ) > 0 ){
				irtmodfit <- TRUE 
				cm1 <- c( cm1 , object$objname )
						}  else { 
				cm1 <- c( cm1 , cl1[1] )
                        }				
		if (irtmodfit){  
			dfr <- t(object$IRT.IC)
			dfr <- as.data.frame(dfr)
			statlist <- object$statlist
			LS <- length(statlist)
			for (ll in 1:LS){
                dfr[ , colnames(statlist)[ll] ] <- 	statlist[,ll]
							}
					} else {
			dfr <- IRT.IC(object)
						}								
						
		LO <- length(object_list)
		#****************************
		#**** loop over remaining objects
		for (vv in 1:LO){
		  irtmodfit <- FALSE
		  if ( length( grep("IRT.modelfit" , class(object_list[[vv]]) ) ) > 0 ){
				irtmodfit <- TRUE 
				obj_vv <- object_list[[vv]]	
				cm1 <- c( cm1 , obj_vv$objname )				
						}	 else {
				cm1 <- c( cm1 , cl1[vv+1] )
							}
			if ( irtmodfit ){                	
				dfr1 <- t(obj_vv$IRT.IC)
				dfr1 <- as.data.frame(dfr1)
				statlist <- obj_vv$statlist
				LS <- length(statlist)
				for (ll in 1:LS){
					dfr1[ , colnames(statlist)[ll] ] <- 	statlist[,ll]
								}
						} else {
				dfr1 <- IRT.IC( object_list[[vv]] )
								}
   		   if ( ! irtmodfit ){  dfr <- rbind( dfr , dfr1 ) }
		   if (irtmodfit ){ 
					dfr <- plyr::rbind.fill( as.data.frame(dfr) , as.data.frame(dfr1 ) ) 
							}
					}
		rownames(dfr) <- NULL

		
		dfr <- data.frame( "Model" = cm1 , dfr )
		IC <- dfr
		# rownames(IC) <- paste(dfr$Model)
		rownames(IC) <- NULL
		
		res <- list("IC" = IC )
		#*************************************
		# collect all likelihood ratio tests
		dfr <- NULL
		M0 <- nrow(IC)		
		for (ii in 1:(M0-1) ){
			for (jj in (ii+1):M0){
				# ii <- 1
				# jj <- 2				
				ii0 <- ii
				jj0 <- jj				
				if ( IC[ii,"Npars" ] > IC[jj,"Npars"] ){ 
						tt <- ii
						ii0 <- jj
						jj0 <- tt
							}				
				if ( IC[ii,"Npars" ] != IC[jj,"Npars"] ){ 							
					dfr1 <- data.frame( "Model1" = IC[ii0,"Model"] , "Model2" = IC[jj0,"Model"]   )
					dfr1$Chi2 <- IC[ii0,"Deviance"] - IC[jj0,"Deviance"]
					dfr1$df <- - ( IC[ii0,"Npars"] - IC[jj0,"Npars"] )
					dfr1$p <- 1 - stats::pchisq( dfr1$Chi2 , df=dfr1$df )
					dfr <- rbind( dfr , dfr1 )
									}
							}
				}
	    res$LRtest <- dfr				
		class(res) <- "IRT.compareModels"
		return(res)
				}
##########################################################################

################################################
# summary method for IRT.compareModels
summary.IRT.compareModels <- function( object , extended = TRUE , ... ){
    dfr1 <- object$IC
	if ( ! extended ){
		vars <- c( "AICc" , "CAIC" , "maxX2" , "MADQ3" , "MADaQ3" ,
						"SRMR" ) 		
		dfr1 <- dfr1[ , ! (colnames(dfr1) %in% vars ) ]		
					}
	cat("Absolute and relative model fit\n\n")
	for ( vv in 2:(ncol(dfr1) ) ){
			dfr1[,vv] <- round( dfr1[,vv] , 3 )
							}
	print(dfr1, ...)			
	dfr2 <- object$LRtest
	if ( ! is.null(dfr2) ){
		cat("\nLikelihood ratio tests - model comparison \n\n")
		obji <- object$LRtest
		for (vv in seq(3,ncol(obji) ) ){
				obji[,vv] <- round( obji[,vv] , 4 )
							}
		print( obji)
	
						}

					}