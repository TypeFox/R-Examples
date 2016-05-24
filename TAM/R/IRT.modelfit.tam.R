
###########################################################
# general model fit function for TAM objects
IRT.modelfit.TAM <- function( object , mod ){
    res <- tam.modelfit( object )
	res$IRT.IC <- IRT.IC(object)
	res$objname <- mod
	class(res) <- paste0("IRT.modelfit." , class(object) )
	return(res)
		}
###########################################################

###########################################################
# IRT.modelfit for objects of class tam.mml
IRT.modelfit.tam.mml <- function( object , ... ){
    cl <- paste(match.call())[2]	
    res <- IRT.modelfit.TAM( object , mod=cl )
	return(res)	
						}
IRT.modelfit.tam.mml.3pl <- IRT.modelfit.tam.mml						
IRT.modelfit.tamaan <- IRT.modelfit.tam.mml	
#####################################################

############################################################
# summary
summary.IRT.modelfit.TAM.helper <- function( object , ... ){	
	cat("Test of Global Model Fit\n")
	obji <- round( object$modelfit.test , 5 )
	print(obji)
	cat("\nFit Statistics\n")
	obji <- object$statlist
	for (vv in seq(1,ncol(obji))){	obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)		
		}
#################################################################	

summary.IRT.modelfit.tam.mml <- summary.IRT.modelfit.TAM.helper
summary.IRT.modelfit.tam.mml.3pl <- summary.IRT.modelfit.TAM.helper
summary.IRT.modelfit.tamaan <- summary.IRT.modelfit.TAM.helper