#*******************************************************
# Summary for tam object                               *
summary.tam.mml <- summary.tam.2pl <- 
	summary.tam.mfr <- summary.tam <- summary.tam.latreg <- 
function( object , file = NULL , ...){


	osink( file = file , prefix = "__SUMMARY.Rout" )
						
    # object      ... object from tam.mml                #
#	object <- tamobject
   latreg <- FALSE
   if ( class(object) == "tam.latreg" ){
		latreg <- TRUE 
		object$irtmodel <- "tam.latreg"
						}

	cat("------------------------------------------------------------\n")
    d1 <- utils::packageDescription("TAM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , 
			sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$time[2] ) , "\n" )
	cat("Computation time:" , print(object$time[2] - object$time[1]), "\n")
	cat( Rsessinfo() , "\n\n")			
    cat("Multidimensional Item Response Model in TAM \n\n")
	irtmodel <- object$irtmodel
	
	cat("IRT Model" , irtmodel ) # , "\n")
	
	# print Call
    print_CALL(object$CALL)	
	
	cat("------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
	
	ctr <- object$control
	if (ctr$snodes==0){
		cat("Numeric integration with" , dim(object$theta)[1] , "integration points\n")
					}
	if (ctr$snodes>0){
		if (ctr$QMC){
			cat("Quasi Monte Carlo integration with" , dim(object$theta)[1] , "integration points\n")
						}
		if (! ctr$QMC){
			cat("Monte Carlo integration with" , dim(object$theta)[1] , "integration points\n")
						}						
					}					
    cat( "\nDeviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$nstud , "\n" )    
    cat( "Number of persons used = " , object$ic$n , "\n" )  
	
	if( ! is.null( object$formulaA)  ){	
	    cat( "Number of generalized items = " , object$nitems , "\n" )    
	    cat( "Number of items = " , ncol(object$resp_orig) , "\n" )  		
								} else {
	    cat( "Number of items = " , object$nitems , "\n" )    				
						}
		
    cat( "Number of estimated parameters = " , object$ic$Npars , "\n" )    
	if (! latreg ){
		cat( "    Item threshold parameters  = " , object$ic$Nparsxsi , "\n" )    
		cat( "    Item slope parameters      = " , object$ic$NparsB , "\n" )    
				  }
    cat( "    Regression parameters      = " , object$ic$Nparsbeta , "\n" )    	
    cat( "    (Co)Variance parameters    = " , object$ic$Nparscov , "\n\n" )    		
	
    cat( "AIC  = " , round( object$ic$AIC , 0 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 0 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 0 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "aBIC  = " , round( object$ic$aBIC , 0 ) , " | penalty =" , round( object$ic$aBIC - object$ic$deviance ,2 ) , 
			"   | aBIC = -2*LL + log((n-2)/24)*p  (adjusted BIC) \n" ) 
    cat( "CAIC = " , round( object$ic$CAIC , 0 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

	cat("------------------------------------------------------------\n")
		cat("EAP Reliability\n")
		obji <- round( object$EAP.rel , 3 )
		print( obji )		
	cat("------------------------------------------------------------\n")
		cat("Covariances and Variances\n")
		if ( object$G >1){
			a1 <- aggregate( object$variance , list( object$group ) , mean )
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
		
	if ( ! latreg ){		
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
			xsi99 <- sum( object$xsi == 99 )
			if ( xsi99 > 0 ){
			  cat("\nSome item xsi parameters are not estimable ")
			  cat(" which is indicated by values of 99\n\n")	
							}
			if ( object$PSF ){
        cat("\nA pseudo facet 'psf' with zero effects with all zero effects\n")
		cat("was created because of non-unique person-facet combinations.\n\n") 							}
							
							
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
		#*******************
		# output efa
		if ( object$irtmodel %in% c("efa") ){
		cat("------------------------------------------------------------\n")	
		cat("\nStandardized Factor Loadings Oblimin Rotation\n")		
		print(object$efa.oblimin)
						}						
		#*******************
		# output bifactor models	
		if ( object$irtmodel %in% c("bifactor1" , "bifactor2","efa") ){
		cat("------------------------------------------------------------\n")	
			if (irtmodel=="efa"){
					cat("\nStandardized Factor Loadings (Schmid Leimann transformation)\n")		
					obji <- object$B.SL
							} else {
					cat("\nStandardized Factor Loadings (Bifactor Model)\n")		
					obji <- object$B.stand				
								}
			meas <- object$meas
			for (vv in seq(1,ncol(obji) ) ){ obji[,vv] <- round( obji[,vv] , 3) }
			print(obji)
			cat("\nDimensionality/Reliability Statistics\n\n")	
			cat("ECV (Omega Asymptotical)=" , round( meas["ECV(omega_a)"] ,3 ) , "\n")
			cat("Omega Total =" , round( meas["omega_t"] ,3 ) , "\n")
			cat("Omega Hierarchical =" , round( meas["omega_h"] ,3 ) , "\n")		
			if (object$maxK==2){		
			cat("Omega Total (GY) =" , round( meas["omega_tot_diff"] ,3 ) , "\n")				
				cat( "  Omega Total GY (Green & Yang, 2009) includes item difficulties\n")
				cat( "  and estimates the reliability of the sum score.\n")	
								}
						}
			}
	#******
	csink(file)
	
			}
#*******************************************************
