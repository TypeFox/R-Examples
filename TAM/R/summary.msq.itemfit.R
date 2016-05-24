


###################################################
# summary for tam.fit
summary.msq.itemfit <- function( object , ... ){

	cat("------------------------------------------------------------\n")
    d1 <- utils::packageDescription("TAM")
	cat( paste( d1$Package , " " , d1$Version , " (" , d1$Date , ")" , 
			sep="") , "\n\n" )	
	cat( "Date of Analysis:" , paste( object$time[2] ) , "\n" )
	cat("Computation time:" , print(object$time[2] - object$time[1]), "\n\n")
	# cat( Rsessinfo() , "\n\n")			

	cat("MSQ item fit statitics (Function 'msq.itemfit')\n\n")

	cat("****************************************************\n")
	cat("\nSummary outfit and infit statistic\n")
    obji <- object$summary_itemfit
	for ( vv in seq(2,ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
					}	
	rownames(obji) <- NULL
	print(obji)
	cat("\n****************************************************\n")
	cat("\nOutfit and infit statistic\n")	
    object <- object$itemfit
	ind <- grep( "fitgroup" , colnames(object) )
	obji <- object
	for ( vv in seq(ind+1,ncol(obji) ) ){
		obji[,vv] <- round( obji[,vv] , 3 )
					}
	print(obji)
	invisible(obji)
		}
###################################################