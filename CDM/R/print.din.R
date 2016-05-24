################################################################################
# print method for objects of class "din"                                      #
################################################################################
print.din <-
function(x, ... ){
	cat("Estimation of Mixed DINA/DINO Model\n\n")

    d1 <- utils::packageDescription("CDM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n" )		
	
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")
	#*** parameters
	cat(paste0("Number of cases = " , x$I , "\n") )
	cat(paste0("Number of items = " , ncol(x$data) , "\n") )
	cat(paste0("Number of skill dimensions = " , ncol(x$q.matrix) , "\n") )	
	cat(paste0("Number of skill classes = " , nrow(x$attribute.patt) , "\n") )		
	cat(paste0("Number of parameters = " , sum(x$Npars) , "\n") )
	cat(paste0("  # item parameters = " , x$Npars$itempars , "\n") )
	cat(paste0("  # skill distribution parameters = " , 
				x$Npars$skillpars , "\n") )
	#*** likelihood
	cat( paste0( "\nLog-Likelihood = " , round( x$loglike ,2 ) , "\n") )
	#*** information criteria
	cat( paste0( "AIC = " , round( x$AIC ,0 ) , "\n") )
	cat( paste0( "BIC = " , round( x$BIC ,0 ) , "\n") )
	invisible(x)	
}

