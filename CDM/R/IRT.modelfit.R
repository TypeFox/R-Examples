
###########################################################
IRT.modelfit <- function (object, ...) {
    UseMethod("IRT.modelfit")
       }
###########################################################
# general model fit function for CDM objects
IRT.modelfit.CDM <- function( object , mod ){
    res <- modelfit.cor.din( dinobj = object)
	res$IRT.IC <- IRT.IC(object)
	res$objname <- mod
	class(res) <- paste0("IRT.modelfit." , class(object) )
	return(res)
		}
###########################################################


###########################################################
# IRT.modelfit for objects of class din, gdina
IRT.modelfit.din <- function( object , ... ){
    cl <- paste(match.call())[2]	
    res <- IRT.modelfit.CDM( object , mod=cl )
	return(res)	
						}
#####################################################
# IRT.modelfit for gdina objects						
IRT.modelfit.gdina <- function( object , ... ){
    cl <- paste(match.call())[2]	
    res <- IRT.modelfit.CDM( object , mod=cl )
	return(res)	
						}						
#############################################################						


############################################################
# summary
summary.IRT.modelfit.helper <- function( object , ... ){	
	cat("Test of Global Model Fit\n")
	obji <- object$modelfit.test
	for (vv in seq(2,ncol(obji))){	obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)
	cat("\nFit Statistics\n")
	obji <- object$modelfit.stat
	for (vv in seq(1,ncol(obji))){	obji[,vv] <- round( obji[,vv] , 3 ) }
	print(obji)		
		}
#################################################################	


# summary.modelfit.cor.din 
summary.IRT.modelfit.din <- summary.IRT.modelfit.helper
summary.IRT.modelfit.gdina <- summary.IRT.modelfit.helper
# summary.IRT.modelfit.gdm <- summary.modelfit.cor.gdm
