
##############################################
# print method for TAM
print_tam <- function( x , ...){
	
	object <- x

    d1 <- utils::packageDescription("TAM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , 
			"" )		

	# print Call
    print_CALL(object$CALL)	

    cat("Multidimensional Item Response Model in TAM \n")
	# irtmodel <- object$irtmodel	
	# cat("IRT Model" , irtmodel  , "\n")	
	
    cat( "\nDeviance = " , round( object$deviance , 2 ) , " | " )
    cat( "Log Likelihood = " , round( -object$deviance/2 , 2 ) , "\n" )	
    # cat( "Number of persons = " , object$nstud , "\n" )    
    cat( "Number of persons used = " , object$ic$n , "\n" )  
    cat( "Number of estimated parameters = " , object$ic$Npars , "\n" )  	

    cat( "AIC  = " , round( object$ic$AIC , 0 ) , "\n" )    
    cat( "BIC  = " , round( object$ic$BIC , 0 ) , "\n" ) 	
	
		}
#########################################################################		

print.tam.mml <- print_tam
# print.tam.mml.2pl <- print_tam
print.tam.mml.3pl <- print_tam
# print.tam.mml.mfr <- print_tam
print.tam.latreg <- print_tam
print.tamaan <- print_tam
print.tam <- print_tam

# summary.tam.mml <- summary.tam.2pl <- 
#	summary.tam.mfr <- summary.tam <- summary.tam.latreg <- 