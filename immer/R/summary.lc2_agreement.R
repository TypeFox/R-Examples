

#*******************************************************
# Summary for immer object
summary.lc2_agreement <- function( object , digits=3 ,... ){
	cat("-----------------------------------------------------------------\n")
	d1 <- utils::packageDescription("immer")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	
	cat( object$description , "\n\n")
	
	cat("Call:\n", paste(deparse(object$CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")			
	
	cat( "Date of Analysis:" , "\n" )
	cat( "   Start:" , paste( object$time$start ) , "\n" )	
	cat( "   End  :" , paste( object$time$end ) , "\n" )		
	cat("Computation time:" , print(object$time$end - object$time$start), "\n\n")    

	# cat("\n")
	
	cat("-----------------------------------------------------------------\n")

	cat( "Deviance = " , round( -2*object$loglike , 2 ) , " | " )
    cat( "Log Likelihood = " , round( object$loglike , 2 ) , "\n" )	

	cat( "Number of persons (sum of weights) =" , object$nobs , "\n" )
    cat( "Number of estimated parameters = " , object$ic$np , "\n\n" )    

		cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
				round( object$ic$AIC - object$ic$dev ,2 ) , 
				"   | AIC = -2*LL + 2*p  \n" )   
		cat( "AIC3 = " , round( object$ic$AIC3 , 2 ) , " | penalty =" , 
					round( object$ic$AIC3 - object$ic$dev ,2 ) , 
					"   | AIC3 = -2*LL + 3*p  \n" )  							
		cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , 
				round( object$ic$AICc - object$ic$dev ,2 ) )
			cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
		cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , 
					round( object$ic$BIC - object$ic$dev ,2 ) , 
				"   | BIC = -2*LL + log(n)*p  \n" )  
		cat( "aBIC = " , round( object$ic$aBIC , 2 ) , " | penalty =" , 
					round( object$ic$aBIC - object$ic$dev ,2 ) , 
				"   | aBIC = -2*LL + log((n-2)/24)*p  (adjusted BIC) \n" ) 
		cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , 
						round( object$ic$CAIC - object$ic$dev ,2 ) )
			cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )  
	
	cat("-----------------------------------------------------------------\n")
	cat("Likelihood Ratio Test of Model Fit \n\n")

	
	cat( "Chi Square = " , round( object$LRT_output$chisquare , 2 ) , " \n" )
    cat( "Degrees of Freedom = " , round( object$LRT_output$df , 2 ) , "\n" )	
	cat( "p-Value = " , round( object$LRT_output$p , 5 ) , "\n\n" )	
	
	cat( "Clogg's NFI = " , round( object$NFI , 3 ) , 
	  "           |  NFI = (L0-L1)/L0 \n\n" )	
	
	cat("-----------------------------------------------------------------\n")
	cat("Model Parameters \n\n")
	
	obji <- paste0( "True agreement: " ,
			"   gamma  = " , round(object$agree_true, digits=digits) )
	cat(obji)
	obji <- paste0( "\nChance agreement: " ,
			" chance = " , round(object$agree_chance, digits=digits) )
	cat(obji)
	cat("\n\n")

	obji <- paste0( "Conditional true agreement: " ,
			" rel = " , round(object$rel_agree, digits=digits) )
	cat(obji)
	cat("\n\n")
	
	obji <- object$parmsummary
	V <- ncol(obji)
	for (vv in 2:V ){ obji[,vv] <- round( obji[,vv] , digits ) }
	rownames(obji) <- NULL
	print(obji)
                }
#*******************************************************
