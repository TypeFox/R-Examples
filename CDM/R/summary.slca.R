

#*******************************************************
# Summary for slca object
summary.slca <- function( object , file = NULL , ... ){
    # object      ... object from rasch.mml                #
	
 	osink( file = file , suffix = paste0( "__SUMMARY.Rout") )

    cat("-----------------------------------------------------------------------------\n")
    d1 <- utils::packageDescription("CDM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )
    
	cat("Call:\n")
	print(object$call)
	cat("\n")
	
	cat( "Date of Analysis:" , paste( object$time$s2 ) , "\n" )
	cat("Computation Time:" , print(object$time$s2 - object$time$s1), "\n\n")
    cat("Structured Latent Class Analysis - Function 'slca' \n")
	modeltype <- object$irtmodel

	cat( "   " , object$N , "Cases, " , object$I , "Items, " , object$G , "Group(s)", "," ,
				object$TP , "Skill classes\n")  
	cat("\n **** Check carefully the number of parameters and identifiability
			of the model.  ***\n")

    if (object$G > 1 ){					
		# group statistics
		cat("\nGroup statistics\n")
		print( object$group.stat )	
					}
    cat("\n-----------------------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
	if ( ! object$converged ){ cat("Maximum number of iterations was reached.\n") }
	cat( "Iteration with minimal deviance =" , object$iter.min , "\n" )
	
	
    cat( "\nDeviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    

    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "  Number of estimated lambda parameters = " , object$ic$itempars , "\n" ) 
    cat( "  Number of estimated distribution parameters = " , object$ic$traitpars , "\n" )    

	
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   


#	cat( "\nEAP Reliability:\n" ) 
#	print( round( t(object$EAP.rel ), 3 ) ) 	
		 cat("-----------------------------------------------------------------------------\n")
		cat("Xlambda Parameters \n")
		obji <- object$Xlambda
		obji <- round( obji , 3)
		print( obji )   
		 cat("-----------------------------------------------------------------------------\n")
		cat("Conditional Item Probabilities \n")
		obji <- object$item
		NV <- ncol(obji)
		for (vv in 3:NV){   obji[,vv] <- round( obji[,vv] , 3) }
		print( obji )   		
		 cat("-----------------------------------------------------------------------------\n")
		cat("Skill Class Parameters \n")
		obji <- object$delta
		obji <- round( obji , 3)
		print( obji ) 		
		 cat("-----------------------------------------------------------------------------\n")
		cat("Skill Class Probabilities \n")
		obji <- object$pi.k
		obji <- round( obji , 4)
		print( obji )   
		
#  cat("\nMean of RMSEA item fit:" , 
#     round( object$mean.rmsea ,3 ) , "\n")		        

   csink( file = file )

			}
#*******************************************************
