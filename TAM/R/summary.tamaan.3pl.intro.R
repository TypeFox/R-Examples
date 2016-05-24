



################################################
# tamaan summary 3pl introduction
summary.tamaan.3pl.intro <- function(object){
		cat("------------------------------------------------------------\n")
		d1 <- utils::packageDescription("TAM")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , 
				sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$time[2] ) , "\n" )
		cat("Computation time:" , print(object$time[2] - object$time[1]), "\n")
		cat( Rsessinfo() , "\n\n")			
		cat("Multidimensional Item Response Model in TAM \n\n")
		irtmodel <- object$irtmodel
		cat("IRT Model" , irtmodel , " (Function 'tam.mml.3pl')\n")
		
		cat("------------------------------------------------------------\n")
		cat( "Number of iterations =" , object$iter , "\n\n" )
		
		ctr <- object$control
		cat("Skill space:" , ifelse(object$skillspace=="normal" ,
				"Normal Distribution", "Discrete Distribution" ) , "\n")	
		if (object$skillspace=="normal"){	
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
								}
		cat( "\nDeviance = " , round( object$deviance , 2 ) , " | " )
		cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
		cat( "Number of persons = " , object$nstud , "\n" )    
		
		if( ! is.null( object$formulaA)  ){	
			cat( "Number of generalized items = " , object$nitems , "\n" )    
			cat( "Number of items = " , ncol(object$resp_orig) , "\n" )  		
									} else {
			cat( "Number of items = " , object$nitems , "\n" )    				
							}
			
		cat( "Number of estimated parameters = " , object$ic$Npars , "\n" )    
		cat( "    Item threshold parameters  = " , object$ic$Nparsxsi , "\n" )    
		cat( "    Item slope parameters      = " , object$ic$NparsB , "\n" )    
		cat( "      Non-active item slopes   = " , 
								object$ic$Ngamma.nonactive , "\n" )    
		cat( "    Item guessing parameters   = " , object$ic$Nguess , "\n" )
		cat( "    Regression parameters      = " , object$ic$Nparsbeta , "\n" )    	
		cat( "    (Co)Variance parameters    = " , object$ic$Nparscov , "\n" )    		
		cat( "    Delta parameters           = " , object$ic$Ndelta , "\n\n" )    	
		
		cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
				round( object$ic$AIC - object$ic$deviance ,2 ) , 
				"   | AIC = -2*LL + 2*p  \n" )    
		cat( "AIC3 = " , round( object$ic$AIC3 , 2 ) , " | penalty =" , 
				round( object$ic$AIC3 - object$ic$deviance ,2 ) , 
				"   | AIC3 = -2*LL + 3*p  \n" )  				
		cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , 
				round( object$ic$AICc - object$ic$deviance ,2 ) )
			cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
		cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , 
					round( object$ic$BIC - object$ic$deviance ,2 ) , 
				"   | BIC = -2*LL + log(n)*p  \n" )  
		cat( "aBIC = " , round( object$ic$aBIC , 2 ) , " | penalty =" , 
					round( object$ic$aBIC - object$ic$deviance ,2 ) , 
				"   | aBIC = -2*LL + log((n-2)/24)*p  (adjusted BIC) \n" ) 
		cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , 
						round( object$ic$CAIC - object$ic$deviance ,2 ) )
			cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   
					}
########################################################