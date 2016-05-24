

#*******************************************************
# Summary for rm.facets object                         *
summary.rm.facets <- function( object , ... ){
	cat("-----------------------------------------------------------------\n")
		d1 <- utils::packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
		cat("Computation time:" , print(object$s2 - object$s1), "\n\n")
    cat("Rater Facet Model with Item/Rater Intercepts and Slopes \n")
	cat("-----------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat( "Number of items   = " , object$ic$VV , "\n" )    		
    cat( "Number of raters  = " , object$ic$RR , "\n" )    	
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "   Number of item parameters   = " , object$ic$np.item , "\n" )    	
    cat( "   Number of rater parameters  = " , object$ic$np.rater , "\n" )    		
    cat( "   Number of distribution parameters  = " , object$ic$np.trait , "\n\n" )    			
			
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

	cat( "Trait Distribution\n" )
	cat( "Mean=" , round( object$mu , 3) , " SD=" , round( object$sigma , 3) ) 
	cat( "\n\nEAP Reliability = ") 
	cat(round( object$EAP.rel,3 ) )
	cat( "\n")
	cat("-----------------------------------------------------------------\n")
	cat("Item Parameters \n")
	obji <- object$item
	rvars <- seq( 2 , ncol(obji ) )
	for (vv in rvars ){ obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)
	cat("-----------------------------------------------------------------\n")
	cat("Rater Parameters \n")
	obji <- object$rater
	rvars <- seq( 2 , ncol(obji ) )
	for (vv in rvars ){ obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)	
                }
#*******************************************************



