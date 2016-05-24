#' Extract a vector of cashflow times in years from a list of instruments
#' 
#' Assumes that LIBOR tenor is in days, with 365 days per year.
#' Assumes that SWAPs are semi-annual
#' Returns a vector of all unique cashflow times in years
#' 
#' @param dfInstruments A dataframe of instuments with at least columns Type and Tenor
#' 
#' 
fCreateTimeVector <- function( dfInstruments ) {
	
	AllTimes <- 0
	
	for ( idx in 1%>%NROW(dfInstruments) ) {
		
		thisInstrument <- dfInstruments[idx,]
		
		if ( "LIBOR" == thisInstrument[["Type"]] ) {
			times <- fGetTimesLibor( thisInstrument )
		} else if ( "SWAP" == thisInstrument[["Type"]] ) {
			times <- fGetTimesSwap(  thisInstrument )
		} else {
			stop( "Unknown instrument Type " %&% thisInstrument[["Type"]] %&% " at line " %&% idx )
		}
	
		AllTimes <- c( AllTimes, times )
		
	}
	
	AllTimes <- AllTimes[-1]	#Remove the redundant zero
	
	UniqueTimes <- sort( unique( AllTimes ) )
	
	return( UniqueTimes )
	
}



#' Extract the payment date of a LIBOR agreement in years
#' 
#' @param dfInstrument A dataframe of instuments with at least columns Type and Tenor
#' 
fGetTimesLibor <- function( dfInstrument ) {
		
	t <- dfInstrument[["Tenor"]]
	return( t )
	
}

#' Extract the payment dates of a Swap agreement in years
#' 
#' @param dfInstrument A dataframe of instuments with at least columns Type and Tenor
#' 
fGetTimesSwap <- function( dfInstrument ) {
	
	freq <- dfInstrument[["Frequency"]]
	t <- sequence( dfInstrument[["Tenor"]] * freq ) / freq
	return( t )
	
}