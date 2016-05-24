#*******************************************************
summary.fuzcluster <- function( object,... ){
    # object      ... object from rasch.mml                #
    cat("---------------------------------------------------------------------------------------------------------- \n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
    cat( "Random starts = " , object$nstarts , "\n" )    
	cat("Computation time:" , print(object$s2 - object$s1), "\n\n")
    cat("---------------------------------------------------------------------------------------------------------- \n")
	cat( "Number of iterations =" , object$iter , 
	      "(Seed =" , object$seed , ")" , "\n" )
    cat( "Deviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
	
    cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , 
			"   | AIC = -2*LL + 2*p  \n" )    
    cat( "AICc = " , round( object$ic$AICc , 2 ) ," | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )
		cat("    | AICc = -2*LL + 2*p + 2*p*(p+1)/(n-p-1)  (bias corrected AIC)\n" )   	
    cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , 
			"   | BIC = -2*LL + log(n)*p  \n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat("   | CAIC = -2*LL + [log(n)+1]*p  (consistent AIC)\n\n" )   

		cat("--------------------------------------------------- \n")
		cat(paste0( "Class probabilities (" , object$K , " Classes)"),
				"\n")
		obji <- object$pi_est
		obji <- round( obji , 4)
		print( obji )                		
		cat("\nMeans\n")
		obji <- object$mu_est
		obji <- round( obji , 4)
		print( obji )   		
		cat("\nStandard deviations\n")
		obji <- object$sd_est
		obji <- round( obji , 4)
		print( obji )   		

			}
#*******************************************************
