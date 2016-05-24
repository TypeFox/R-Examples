



########################################################
summary.tamaan.3pl.lcaprobs <- function(object){
		cat("------------------------------------------------------------\n")
		cat("Item Response Probabilities\n")
		obji <- object$lcaprobs
		obji[,-1] <- round( obji[,-1] , 4 )	
		print(obji)
			}
##########################################################


############################################################
# discrete distributions
summary.tamaan.3pl.discrete.distribution <- function(object){
			# latent class distributions
			cat("------------------------------------------------------------\n")
			cat("Trait distribution parameters delta\n")
			obji <- round( object$delta , 4 )
			colnames(obji) <- paste0("Group" , 1:object$G)
			print( obji )		
			TP <- nrow(obji)			
			cat("------------------------------------------------------------\n")
			cat("Full Trait distribution\n")
			obji <- round( object$pi.k , 4 )			
			colnames(obji) <- paste0("Group" , 1:object$G)
			print( obji )		
			}
############################################################	



###############################################################
# normal skill space
summary.tamaan.normal.skillspace <- function(object){
		cat("------------------------------------------------------------\n")
			cat("EAP Reliability\n")
			obji <- round( object$EAP.rel , 3 )
			print( obji )		
		cat("------------------------------------------------------------\n")
			cat("Covariances and Variances\n")
			if ( object$G >1){
				a1 <- stats::aggregate( object$variance , list( object$group ) , mean )
				object$variance <- a1[,2]
						}
			obji <- round( object$variance , 3 )
			if ( object$G >1){
	#			names(obji) <- paste0("Group" , seq(1,object$G) )
				names(obji) <- paste0("Group" , object$groups )
						}		
			print( obji )
		cat("------------------------------------------------------------\n")
			cat("Correlations and Standard Deviations (in the diagonal)\n")
			if ( object$G >1){
				obji <- sqrt( object$variance )
						} else {
			obji <- stats::cov2cor(object$variance)
			diag(obji) <- sqrt( diag( object$variance) )
						}
			if ( object$G >1){
	# 		names(obji) <- paste0("Group" , seq(1,object$G) )
				names(obji) <- paste0("Group" , object$groups )			
						}		
			obji <- round( obji, 3 )
			print( obji )
		cat("------------------------------------------------------------\n")
		 cat("Regression Coefficients\n")
			obji <- round( object$beta , 5 )
			print( obji )		
				} 
######################################################################				


#########################################################################
# item parameters
summary.tamaan.item.parameters <- function(object){
	cat("------------------------------------------------------------\n")		
		cat("Item Parameters -A*Xsi\n")
#		cat("   Item difficulties -A*Xsi are displayed in 'AXsi_'! \n\n")
		obji <- object$item
		for (vv in seq(2,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
		print(obji)
	# print xsi parameters if 
	if( ! is.null( object$formulaA)  ){
		cat("\nItem Facet Parameters Xsi\n")
#		cat("   Item difficulties -A*Xsi are displayed in 'AXsi_'! \n\n")
		obji <- object$xsi.facets
		for (vv in seq(3,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
		print(obji)
					}				
	if (( object$maxK > 2 ) | ( object$printxsi) ){
		cat("\nItem Parameters Xsi\n")
#		cat("   Item difficulties -A*Xsi are displayed in 'AXsi_'! \n\n")
		obji <- object$xsi
#		obji[,1] <- obji[,-1]
		for (vv in seq(1,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
		print(obji)
			}
		
	cat("------------------------------------------------------------\n")		
	if (( object$maxK > 2 ) | ( object$printxsi) ){
		cat("\nItem Parameters Xsi\n")
#		cat("   Item difficulties -A*Xsi are displayed in 'AXsi_'! \n\n")
		obji <- object$xsi
#		obji[,1] <- obji[,-1]
		for (vv in seq(1,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
		print(obji)
			}
		cat("\nGammaslope Parameters\n")
		obji <- object$gammaslope
		print(round(obji,3)  )   			
		cat("\nGuessing Parameters\n")
		obji <- object$item$guess
		names(obji) <- colnames(object$resp)
		print(round(obji,3)  )   					
			}
###################################################################################					
		


#########################################################################
# item parameters
summary.tamaan.item.parameters.mixture <- function(object){

		cat("------------------------------------------------------------\n")
		cat("Item Parameters\n")
		obji <- object$itempartable_MIXTURE
		for ( vv in 3:(ncol(obji) ) ){
				obji[,vv] <- round( obji[,vv] , 3 )
						}
		print( obji )

		cat("------------------------------------------------------------\n")
		cat("\nGammaslope Parameters\n")
		obji <- object$gammaslope
		print(round(obji,3)  )   								
			}
###################################################################################					
		
		

##############################################res#################
# cluster locations
summary.tamaan.3pl.loclca <- function(object){
		cat("*******************************\n")
		cat("Cluster locations\n")
		obji <- round( object$locs, 3 )
		print( obji )
		#*****************
		# item parameters
#		cat("------------------------------------------------------------\n")		
#		cat("Item Parameters\n")	
#		obj1 <- object$tamaanify$loclca_ITEMS
#		N1 <- nrow(obj1)
#		ipars <- object$gammaslope[ 1:N1 ]
#		ind <- match( paste(obj1$parm) , names(ipars) )
#		obj1$est <- round( ipars , 3 )[ind]
#		obji <- obj1
#		print( obji)		
		# item response probabilities
		summary.tamaan.3pl.lcaprobs(object)
				} 
######################################################################				


###################################################################
# distribution mixture
summary.tamaan.3pl.distr.mixture <- function(object){
		cat("------------------------------------------------------------\n")	
		cat("Class Probabilities\n")
		obji <- round( object$probs_MIXTURE , 3 )
		print( obji )
		cat("------------------------------------------------------------\n")	
		cat("Class Distributions\n\n")
		mom <- object$moments_MIXTURE
		ncl <- length(mom)
		for (cl in 1:ncl){
		   cat("******\nClass" , cl , "\n\n")
		   mom.cl <- mom[[cl]]
		   mom.cl$skewness.trait <- NULL	   
		   for (nn in names(mom.cl)[1:2]){
				mom.cl[[nn]] <- round( mom.cl[[nn]] , 3 )
									}
		   for (nn in names(mom.cl[[3]]) ){
				mom.cl[[3]][[nn]] <- round( mom.cl[[3]][[nn]] , 3 )
									}																		
        	print( mom.cl )		   	 
						}
		   				
			}
#################################################################			