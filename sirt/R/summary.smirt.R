#*******************************************************
# Summary for smirt object
summary.smirt <- function( object , ...){
	cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    cat("  Function 'smirt' \n")
	if (object$irtmodel == "noncomp"){
		cat("   Noncompensatory item response model\n\n")
					}
	if (object$irtmodel == "comp"){
		cat("   Compensatory item response model\n\n")
					}				
	if (object$irtmodel == "partcomp"){
		cat("   Partially compensatory item response model\n\n")
					}
					
#	modeltype <- object$irtmodel
		cat( "   " , object$ic$n , "Cases, " , ncol(object$dat2) , "Items, " , # object$G , "Group(s)", "," ,
				object$D , "Dimension(s)\n")  

    cat("-----------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    

    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "  Number of estimated item parameters = " , object$ic$np.item , "\n" ) 
    cat( "               b parameters = " , object$ic$np.item.b , "\n" ) 	
    cat( "               a parameters = " , object$ic$np.item.a , "\n" ) 		
    cat( "               c parameters = " , object$ic$np.item.c , "\n" ) 
    cat( "               d parameters = " , object$ic$np.item.d , "\n" ) 			
	cat( "            mu.i parameters = " , object$ic$np.item.mu.i , "\n" )
    cat( "  Number of estimated distribution parameters = " , object$ic$np.cov , 
				"\n" )    
    cat( "               Means = " , object$ic$np.cov.mu , "\n" ) 
    cat( "               Covariances = " , object$ic$np.cov.covM , "\n" ) 	
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

    cat("-----------------------------------------------------------------\n")
#	cat("Trait Distribution\n")
#	obji <- object$pi.k
	cat( "\nM Trait:\n" ) 
	print( round( object$mean.trait , 3 )  )
	cat( "\nSD Trait:\n" ) 
	print( round( object$sd.trait , 3 ) ) 
	cat( "\nCorrelations Trait: \n" )
	print( round( object$cor.trait , 3 ) )
	cat( "\nEAP Reliability:\n" ) 
	print( round( t(object$EAP.rel ), 3 ) ) 	
    cat("-----------------------------------------------------------------\n")
	cat("Item Parameter \n")
	obji <- object$item
	obji[,-1] <- round( obji[,-1] , 3)
	print( obji )   
#  cat("\nMean of RMSEA item fit:" , 
#     round( object$mean.rmsea ,3 ) , "\n")		        
			}
#*******************************************************
