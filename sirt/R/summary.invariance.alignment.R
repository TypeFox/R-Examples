#*******************************************************
# Summary for smirt object
summary.invariance.alignment <- function( object , ...){
	cat("-----------------------------------------------------------------\n")
    d1 <- utils::packageDescription("sirt")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
	cat("Computation Time:" , print(object$s2 - object$s1), "\n\n")
    cat("  Function 'invariance.alignment' \n")

    cat("-----------------------------------------------------------------\n")
#	cat( "Number of iterations =" , object$Niter[1] , 
#			"(LAMBDA) | " , object$Niter[2] , "(NU)\n" )
#	cat( "Iteration with minimum function value =" , 
#		object$miniter[1] , 	"(LAMBDA) | " , object$miniter[2] , "(NU)\n" )
    cat( "Optimization Function Value (minimum value) = " , 
		round(object$fopt[1],4) , 	"(LAMBDA) | " , round(object$fopt[2],4) , "(NU)\n" )
#    cat( "Optimization Function Value (last value) = " , 
#				object$fopt.history[object$Niter[1]-1,1] , 	"(LAMBDA) | " , 
#				object$fopt.history[object$Niter[2]-1,2] , "(NU)\n" )
#    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
#    cat( "Number of persons = " , object$ic$n , "\n" )    
    cat("-----------------------------------------------------------------\n")
	cat("Effect Sizes of Approximate Invariance \n")
    print( round( object$es.invariance , 4 ))
    cat("-----------------------------------------------------------------\n")
	cat("Group Means and Standard Deviations \n")
	obji <- object$pars
	obji <- round( obji , 3)
	print( obji )   
    cat("-----------------------------------------------------------------\n")
	cat("Summary Aligned Item Parameters \n")
	obji <- object$itempars.aligned
	obji <- round( obji , 3)
	print( obji )   
    cat("-----------------------------------------------------------------\n")
	cat("Summary Absolute Residuals Loadings LAMBDA \n")
	obji <- as.vector( abs( object$lambda.resid  ))
	obji <- round( summary(obji) , 4)
	print( obji )   
    cat("-----------------------------------------------------------------\n")
	cat("Summary Absolute Residuals Intercepts NU \n")
	obji <- as.vector( abs( object$nu.resid  ))
	obji <- round( summary(obji) , 4)
	print( obji )   

			}
#*******************************************************
