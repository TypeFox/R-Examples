#*******************************************************
# Summary for rasch.copula object                         *
summary.rasch.copula3 <- function( object , ... ){
    # object      ... object from rasch.copula                #
        cat("-----------------------------------------------------------------\n")
		d1 <- utils::packageDescription("sirt")
#		cat( paste( d1$Package , d1$Version ,d1$Date ) , "\n" )
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	    cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")		
        cat("Marginal Maximum Likelihood Estimation \n")
        cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
		cat("Function rasch.copula3\n")
		cat("alpha1=",round(object$alpha1,3)," alpha2=" , round(object$alpha2,3) , " \n")
		cat("-----------------------------------------------------------------\n")
		cat( "Deviance = " , round( object$deviance , 2 ) , "\n" )
		cat( "Number of persons = " , object$ic$n , " (" , nrow(object$pattern) , " Response Patterns)\n" )    
		cat( "Number of estimated parameters = " , object$ic$np , "\n" )   
		cat( "Number of iterations = " , object$iter , "\n" )   	
		cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
					round( object$ic$AIC - object$ic$deviance ,2 ) , "\n" )    
		cat( "AICc = " , round( object$ic$AICc , 2 ) , " | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )    
		cat(" (bias corrected AIC)\n" )   	
		cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , "\n" )  
		cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat( " (consistent AIC) \n\n" )         
		#******		
		cat( "Trait Distribution (" , length(object$theta.k) , " Knots )\n" )
		cat("\nMean Vector\n")
		print( round( object$mu ,3 ))		
		cat("\nCorrelation Matrix\n")
		print( round( object$sigma ,4 ))
		
		cat(paste("\nEAP Reliability:\n "))
		print( round( object$EAP.Rel,3))
		
		
		
		cat("-----------------------------------------------------------------\n")
		cat("Item Parameter \n")
		.pr( object$item , 3 )   
		cat("\nDependency parameters\n")
		.pr(object$summary.delta , digits = 3)
                }
#*******************************************************




#*******************************************************
# Summary for rasch.copula object                         *
summary.rasch.copula2 <- function( object , ... ){
    # object      ... object from rasch.copula                #
        cat("-----------------------------------------------------------------\n")
		d1 <- utils::packageDescription("sirt")
#		cat( paste( d1$Package , d1$Version ,d1$Date ) , "\n" )
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	    cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")		
        cat("Marginal Maximum Likelihood Estimation \n")
        cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
		cat("Function rasch.copula2\n")
		cat("alpha1=",round(object$alpha1,3)," alpha2=" , round(object$alpha2,3) , " \n")
		cat("-----------------------------------------------------------------\n")
		cat( "Deviance = " , round( object$deviance , 2 ) , "\n" )
		cat( "Number of persons = " , object$ic$n , " (" , nrow(object$pattern) , " Response Patterns)\n" )    
		cat( "Number of estimated parameters = " , object$ic$np , "\n" )   
		cat( "Number of iterations = " , object$iter , "\n" )   	
		cat( "AIC  = " , round( object$ic$AIC , 2 ) , " | penalty =" , 
					round( object$ic$AIC - object$ic$deviance ,2 ) , "\n" )    
		cat( "AICc = " , round( object$ic$AICc , 2 ) , " | penalty =" , round( object$ic$AICc - object$ic$deviance ,2 ) )    
		cat(" (bias corrected AIC)\n" )   	
		cat( "BIC  = " , round( object$ic$BIC , 2 ) , " | penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , "\n" )  
		cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," | penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) )
		cat( " (consistent AIC) \n\n" )         
		#******		
		cat( "Trait Distribution (" , length(object$theta.k) , " Knots )\n" , 
				  "Mean=" , 0 , " SD=" , 1 , "\n") 
		cat(paste("\nEAP Reliability:" , round( object$EAP.Rel,3)),"\n\n")			  
		cat("-----------------------------------------------------------------\n")
		cat("Item Parameter \n")
		.pr( object$item , 3 )   
		cat("\nDependency parameters\n")
		.pr(object$summary.delta , digits = 3)
                }
#*******************************************************




#***************************************************
# This is an auxiliary function which helps for 
#   printing some progress
.pr <-function( object,digits){	
			ow <- options()$warn
			if ( length(dim(object)) == 2 ){
					options( warn = -1 )
					if ( nrow(object) >= 1 ){ 
						g1a <- apply( object , 2 , as.numeric )
								} else { g1a <- object }				
					g1a <- matrix(g1a , nrow= nrow(object) , ncol= ncol(object))				
					colnames(g1a) <- colnames(object)
					g1 <- colMeans( g1a )
					g1 <- which( ! is.na( g1 ) )
					options( warn = ow )
					object1 <- object
					object1[ , g1 ] <- round( object1[ , g1 ] , digits )
					# print( object1)					
					print( object1 )  } 
					else 
					{ print( round( object , digits ) ) }
				}
#***************************************************



#*******************************************************
# Summary for rasch.copula object                         *
##NS S3method(summary,rasch.copula)
summary.rasch.copula <- function( object , ...){
    # object      ... object from rasch.copula                #
        cat("---------------------------------------------------------------------------------------------------------- \n")
        cat("Marginal Maximum Likelihood Estimation \n")
        cat(paste( "Raschtype Copula Model with generalized logistic link function: Estimation of alpha1 and alpha2 \n") )
		cat("alpha1=",round(object$alpha1,3)," alpha2=" , round(object$alpha2,3) , " \n")
		cat("Function rasch.copula\n")		
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat( "Deviance = " , round( object$deviance , 2 ) , "\n" )
    cat( "Number of persons = " , object$ic$n , " (" , nrow(object$pattern) , " Response Patterns)\n" )    
    cat( "Number of estimated parameters = " , object$ic$np , "\n" )    
    cat( "AIC = " , round( object$ic$AIC , 2 ) , " ; penalty =" , round( object$ic$AIC - object$ic$deviance ,2 ) , "\n" )    
    cat( "BIC = " , round( object$ic$BIC , 2 ) , " ; penalty =" , round( object$ic$BIC - object$ic$deviance ,2 ) , "\n" )  
    cat( "CAIC = " , round( object$ic$CAIC , 2 ) ," ; penalty =" , round( object$ic$CAIC - object$ic$deviance ,2 ) , "\n\n" )         
    cat( "Trait Distribution (" , length(object$theta.k) , " Knots )\n" , 
              "Mean=" , 0 , " SD=" , 1 , "\n") 
    cat("---------------------------------------------------------------------------------------------------------- \n")
    cat("Item Parameter \n")
    .pr( object$item , 3 )                
                }
#*******************************************************



