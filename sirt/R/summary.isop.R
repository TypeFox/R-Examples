

#*******************************************************
# Summary for isop object                         *
summary.isop <- function( object , ... ){
	cat("-----------------------------------------------------------------\n")
		d1 <- utils::packageDescription("sirt")
		cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , sep="") , "\n\n" )	
		cat( "Date of Analysis:" , paste( object$s2 ) , "\n" )
    cat("ISOP and ADISOP Model \n")
	cat("-----------------------------------------------------------------\n")
    cat( "Number of persons = " , dim(object$dat)[1] , "\n" )    
    cat( "Number of items = " , dim(object$dat)[2] , "\n\n" )    
	
    cat( "Number of person score groups = " , dim(object$freq.correct)[1] , "\n" )    
    cat( "Number of item score groups = " , dim(object$freq.correct)[2] , "\n" )    

	cat( "\nLog-Likelihood Comparison\n") 
	obji <- object$ll
	obji[,2] <- round( obji[,2] , 3 )
	obji[,3] <- round( obji[,3] , 5 )	
	if ( object$model == "np.dich" ){
		obji$fit <- round( c(NA,object$fit ) , 4 )
					}
	print(obji)	
	
	cat("*****************************************************************\n")
	cat("Item Statistics and Scoring\n")
	obji <- object$item
	for (vv in 2:( ncol(obji))){
		obji[,vv] <- round( obji[,vv] , 4 )
					}
	print(obji)
	
                }
#*******************************************************



