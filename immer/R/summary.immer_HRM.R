

#*******************************************************
# Summary for immer object
summary.immer_HRM <- function( object , digits=3 ,... ){
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

	cat( "Number of iterations =" , object$iter , "\n" )
	cat( "Number of burnin iterations =" , object$burnin , "\n" )			
	cat( "Number of saved iterations =" , object$N.save , "\n\n" )			
	
	cat("-----------------------------------------------------------------\n")
    # cat( "Dbar = " , round( object$ic$Dbar , 2 ) , "\n" )# , " | " )
    # cat( "Dhat = " , round( object$ic$Dhat , 2 ) , "\n" )# , " | " )	
    # cat( "pD   = " , round( object$ic$pD , 2 ) ,  " | pD = Dbar - Dhat \n" )		
    # cat( "DIC  = " , round( object$ic$DIC , 2 ) ,  " | DIC = Dhat + pD \n\n" )				
#    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	    

	cat( "Deviance = " , round( object$ic$dev , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$ic$dev/2 , 2 ) , "\n" )	
	cat( "Number of person-rater-interactions = " , object$ic$ND , "\n" )    
	cat( "Number of persons   = " , object$ic$N , "\n" )	
	cat( "Number of items   = " , object$ic$I , "\n" )	
	cat( "Number of raters   = " , object$ic$R , "\n\n" )	
	
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "                      # mu     = " , object$ic$Npars["mu"] , "\n" )    
    cat( "                      # sigma  = " , object$ic$Npars["sigma"] , "\n" )    
    cat( "                      # a      = " , object$ic$Npars["a"] , "\n" )    	
	cat( "                      # b      = " , object$ic$Npars["b"] , "\n" )
    cat( "                      # phi    = " , object$ic$Npars["phi"] , "\n" )    		
	cat( "                      # psi    = " , object$ic$Npars["psi"] , "\n\n" )    		

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
	cat( "Trait Distribution\n" )
	cat( "Mean=" , round( object$mu , 3) , " SD=" , round( object$sigma , 3) , "\n\n") 
	cat( "EAP Reliability = ") 
	cat(round( object$EAP.rel,3 ) )
	cat( "\n")	
	
	cat("-----------------------------------------------------------------\n")
	cat("Item Parameters \n")
	obji <- object$item
	V <- ncol(obji)
	for (vv in 2:V ){ obji[,vv] <- round( obji[,vv] , digits ) }
	rownames(obji) <- NULL
	print(obji)
	cat("-----------------------------------------------------------------\n")
	cat("Rater Parameters \n")
	obji <- object$rater_pars
	V <- ncol(obji)
	for (vv in 4:V ){ obji[,vv] <- round( obji[,vv] , digits ) }
	rownames(obji) <- NULL
	print(obji)
	cat("-----------------------------------------------------------------\n")
	cat("MCMC Diagnostics \n")
	summ <- object$summary.mcmcobj
	obji <- summ[ , c("parameter" , "Mean" , "MAP" , "skewness" ,"SD" ) ] 
	obji$MCMC_SE <- summ$Time.series.SE
	vars <- c("Rhat","PercSERatio","effSize" ) # ,"sampSize")
	for (vv in vars ){
	  obji[,vv] <- summ[,vv]
					}
	V <- ncol(obji)
	for (vv in 2:V ){ obji[,vv] <- round( obji[,vv] , digits ) }		
	rownames(obji) <- NULL
	print(obji)
                }
#*******************************************************
