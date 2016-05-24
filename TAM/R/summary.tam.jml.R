#*******************************************************
# Summary for tam object                 *
summary.tam.jml <- function( object , file = NULL , ...){
	if ( ! is.null( file ) ){
		sink( paste0( file , "__SUMMARY.Rout") , split=TRUE )
						}
    # object      ... object from tam.mml                #
#	object <- tamobject
	cat("------------------------------------------------------------\n")
    d1 <- utils::packageDescription("TAM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , 
			sep="") , "\n\n" )	
	cat( "Start of Analysis:" , paste( object$time[1] ) , "\n" )
#	cat( "Date of Analysis:" , paste( object$time["s2"] ) , "\n" )
	cat( "End of Analysis:" , paste( object$time[2] ) , "\n" )
# print( object$time[2] - object$time[1] )
	cat("Computation time:" , print( object$time[2] - object$time[1] ) , "\n")
	cat( Rsessinfo() , "\n\n")			
    cat("Joint Maximum Likelihood Estimation in TAM \n\n")
	irtmodel <- object$irtmodel
	cat("IRT Model" , irtmodel , "\n")
	
	cat("------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
					
    cat( "\nDeviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$nstud , "\n" )    
	
	if( ! is.null( object$formulaA)  ){	
	    cat( "Number of generalized items = " , object$nitems , "\n" )    
	    cat( "Number of items = " , ncol(object$resp_orig) , "\n" )  		
								} else {
	    cat( "Number of items = " , object$nitems , "\n" )    				
						}
		
#    cat( "Number of estimated parameters = " , object$ic$Npars , "\n" )    
#    cat( "    Item threshold parameters  = " , object$ic$Nparsxsi , "\n" )    
#    cat( "    Item slope parameters      = " , object$ic$NparsB , "\n" )    
#    cat( "    Regression parameters      = " , object$ic$Nparsbeta , "\n" )    	
#    cat( "    (Co)Variance parameters    = " , object$ic$Nparscov , "\n\n" )    		
	
#    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
#			"   | AIC = -2*LL + 2*p  \n" )    
#    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
#		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
#    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
#			"   | BIC = -2*LL + log(n)*p  \n" )  
#    cat( "aBIC  = " , round( object$ic$aBIC , 2 ) , " | penalty =" , round( object$ic$aBIC - object$ic$deviance ,2 ) , 
#			"   | aBIC = -2*LL + log((n-2)/24)*p  (adjusted BIC) \n" ) 
#    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
#		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

#	cat("------------------------------------------------------------\n")
#		cat("EAP Reliability\n")
#		obji <- round( object$EAP.rel , 3 )
#		print( obji )		
#	cat("------------------------------------------------------------\n")
#		cat("Covariances and Variances\n")
#		if ( object$G >1){
#			a1 <- aggregate( object$variance , list( object$group ) , mean )
#			object$variance <- a1[,2]
#					}
#		obji <- round( object$variance , 3 )
#		if ( object$G >1){
#			names(obji) <- paste0("Group" , seq(1,object$G) )
#			names(obji) <- paste0("Group" , object$groups )
#					}		
#		print( obji )
#	cat("------------------------------------------------------------\n")
#		cat("Correlations and Standard Deviations (in the diagonal)\n")
#		if ( object$G >1){
#			obji <- sqrt( object$variance )
#					} else {
#		obji <- cov2cor(object$variance)
#		diag(obji) <- sqrt( diag( object$variance) )
#					}
#		if ( object$G >1){
# 		names(obji) <- paste0("Group" , seq(1,object$G) )
#			names(obji) <- paste0("Group" , object$groups )			
#					}		
#		obji <- round( obji, 3 )
#		print( obji )
#	cat("------------------------------------------------------------\n")
#	 cat("Regression Coefficients\n")
#        obji <- round( object$beta , 5 )
#		print( obji )		
#	cat("------------------------------------------------------------\n")		
#		cat("Item Parameters -A*Xsi\n")
#		cat("   Item difficulties -A*Xsi are displayed in 'AXsi_'! \n\n")
#		obji <- object$item
#		for (vv in seq(2,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
#		print(obji)
#	# print xsi parameters if 
#	if( ! is.null( object$formulaA)  ){
#		cat("\nItem Facet Parameters Xsi\n")
#		cat("   Item difficulties -A*Xsi are displayed in 'AXsi_'! \n\n")
#		obji <- object$xsi.facets
#		for (vv in seq(3,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
#		print(obji)
#					}				
#	if (( object$maxK > 2 ) | ( object$printxsi) ){
		cat("\nItem Parameters Xsi\n")
#		cat("   Item difficulties -A*Xsi are displayed in 'AXsi_'! \n\n")
		obji <- object$item
#		obji[,1] <- obji[,-1]
		for (vv in seq(2,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
		print(obji)
#			}
	#******
	if ( ! is.null( file ) ){
		sink(  )
						}
			}
#*******************************************************
