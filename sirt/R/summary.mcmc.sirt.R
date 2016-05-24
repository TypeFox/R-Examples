

#*******************************************************
# Summary for rm.facets object                         *
summary.mcmc.sirt <- function( object , digits=3 ,... ){
	cat("-----------------------------------------------------------------\n")
	d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , "\n" )
	cat( "   Start:" , paste( object$time$start ) , "\n" )	
	cat( "   End  :" , paste( object$time$end ) , "\n" )		
	cat("Computation time:" , print(object$time$end - object$time$start), "\n\n")    

	cat( "Number of iterations =" , object$iter , "\n" )
	cat( "Number of burnin iterations =" , object$burnin , "\n\n" )			
	cat( object$description , "\n")
	cat("-----------------------------------------------------------------\n")
    cat( "Dbar = " , round( object$ic$Dbar , 2 ) , "\n" )# , " | " )
    cat( "Dhat = " , round( object$ic$Dhat , 2 ) , "\n" )# , " | " )	
    cat( "pD   = " , round( object$ic$pD , 2 ) ,  " | pD = Dbar - Dhat \n" )		
    cat( "DIC  = " , round( object$ic$DIC , 2 ) ,  " | DIC = Dhat + pD \n\n" )				
#    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	    
	cat( "Number of persons = " , nrow(object$dat) , "\n" )    
	if ( object$model == "2pno.ml"){
		cat( "Number of groups  = " , object$ic$G , "\n")
		cat( "  Group sizes: M =" , round(object$ic$M.n,3) , 
			" | SD =" , round(object$ic$SD.n,3) , "\n")
				}
    cat( "Number of items   = " , ncol(object$dat) , "\n" )
	
	
#    cat( "Number of estimated parameters = " , object$np , "\n" )    	
#    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
#			"   | AIC = -2*LL + 2*p  \n" )    
#    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
#		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
#    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
#			"   | BIC = -2*LL + log(n)*p  \n" )  
#    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
#		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

#	cat( "Trait Distribution\n" , 
#			  "Mean=" , 0 , " SD=" , round( object$sigma , 3) ) 
#			  "Mean=" , 0 , " SD=" , 1 ) 

	cat( "\nEAP Reliability = ") 
	cat(round( object$EAP.rel,3 ) )
	cat( "\n")
	cat("-----------------------------------------------------------------\n")
	cat("Item Parameters \n")
	obji <- object$summary.mcmcobj
	obji$PercSEratio <- obji$PercSERatio
	vars <- c("parameter" , "Mean" , "SD" , "MAP" , "Rhat" , "PercSEratio" ,
#			"Q2.5" , "Q5" , "Q95" , "Q97.5" )
			"Q5" , "Q95"  )
	obji$PercSEratio <- round( obji$PercSEratio	, 1 )
	obji$Rhat <- round( obji$Rhat	, 2 )
	obji <- obji[ , vars ]
	rvars <- seq( 2 , ncol(obji ) )
	for (vv in rvars ){ obji[,vv] <- round( obji[,vv] , digits ) }
	rownames(obji) <- NULL
	print(obji)
                }
#*******************************************************



