#*******************************************************
# Summary for gom object
summary.gom <- function( object ,...){
	cat("----------------------------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    cat("  Function 'gom.em' \n")
	if (object$model == "GOM"){
		cat("   Discrete Grade of Membership Model\n\n")
					}
	if (object$model == "GOMRasch"){
		cat("   Rasch Grade of Membership Model\n\n")
					}					
	modeltype <- object$irtmodel
		cat( "   " , object$ic$n , "Cases, " , ncol(object$dat2) , "Items, " , 
				object$K , "Classes", "," ,
				object$TP , "Discrete Integration Points\n")  

	cat("----------------------------------------------------------------------------------\n")
	cat( "Number of iterations =" , object$iter , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of persons = " , object$ic$n , "\n" )    

    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "  Number of estimated item parameters = " , object$ic$np.item , "\n" ) 
    cat( "  Number of estimated distribution parameters = " , object$ic$np.trait , 
				"\n" )    
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
				round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , 
				round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , 
			round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , 
			round( object$ic$CAIC - object$ic$deviance ,2 ) )
	cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   
	#*******************************************
    if (object$model=="GOMRasch"){
		cat("----------------------------------------------------------------------------------\n")
		cat("Trait Distribution (Location, Variability)\n")		
		cat( " Means: " , round( object$mu , 3 ) , "\n")				
		cat( " Standard deviations: " , round( sqrt(diag(object$Sigma)) , 3 ) , "\n")
		c1 <- cov2cor(object$Sigma)
		cat( " Correlation " , round( c1[lower.tri(c1)] , 3 ) , "\n")		
		cat("EAP Reliability = " , round(object$EAP.rel,3) , "\n")			
			}
	#*******************************************			
    if ( ! ( object$plmat) ){	
		cat("----------------------------------------------------------------------------------\n")
		cat("Membership Function Descriptives \n")
		obji <- object$classdesc
		obji <- round( obji , 3)
		print( obji )   
				}
	cat("----------------------------------------------------------------------------------\n")
	cat("Item Parameters \n")
	obji <- object$item
	obji[,-1] <- round( obji[,-1] , 3)
	print( obji )   
#  cat("\nMean of RMSEA item fit:" , 
#     round( object$mean.rmsea ,3 ) , "\n")		        
			}
#*******************************************************
