
###################################################
# summary for objects of class msq.itemfitWLE
summary.msq.itemfitWLE <- function( object , ... ){

	cat("------------------------------------------------------------\n")
    d1 <- utils::packageDescription("TAM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , 
			sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$time[2] ) , "\n" )
	cat("Computation time:" , print(object$time[2] - object$time[1]), "\n\n")
	# cat( Rsessinfo() , "\n\n")			

	cat("MSQ item fit statitics (Function 'msq.itemfitWLE')\n\n")

	cat("****************************************************\n")
	cat("\nSummary outfit and infit statistic\n")
	
	if ( is.null(object$fitindices) ){
		object1 <- object$fit_data_summary 
				} else {
		object1 <- object$fit_parm_summary 
					}	
	
    obji <- object1
	for ( vv in seq(2,ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
					}	
	rownames(obji) <- NULL
	print(obji)
	
	cat("\n****************************************************\n")
	cat("\nOutfit and infit statistic\n")		
	
	
	if ( is.null(object$fitindices) ){
		object <- object$fit_data 
				} else {
		object <- object$fit_parm 
					}
	obji <- object
	for ( vv in seq(2,ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
					}
    rownames(obji) <- NULL					
	print(obji)
		}
###################################################