#' Returns the matrix of cashflows for the list of instruments
#' 
#' @param dfInstruments A set of market instruments as a dataframe with columns Type, Tenor, Frequency and Rate with Type in (LIBOR, SWAP), Tenor the instrument maturity in years and rate the rate per annum
#' 
#' 
fCreateCashflowMatrix <- function( dfInstruments ) {

	Times <- fCreateTimeVector( dfInstruments )
	
	Cashflows <- matrix( 0, nrow= NROW( dfInstruments ), ncol= length( Times ) )
	
	for ( idx in 1%>%NROW(dfInstruments) ) {
		
		thisInstrument <- dfInstruments[idx,]
		
		if ( "LIBOR" == thisInstrument[["Type"]] ) {
			cashflowSchedule <- fGetCashflowsLibor( thisInstrument )
		} else if ( "SWAP" == thisInstrument[["Type"]] ) {
			cashflowSchedule <- fGetCashflowsSwap( thisInstrument )
		} else {
			stop( "Unknown instrument Type " %&% thisInstrument[["Type"]] %&% " at line " %&% idx )
		}
		
		Cashflows[idx, match( cashflowSchedule$times, Times ) ] <- cashflowSchedule$cashflows
	}
	
	colnames( Cashflows ) <- Times
	return( Cashflows )
	
}


#' Gets the cashflow schedule for a LIBOR agreement
#' 
#' @param dfInstrument A set of market instruments as a dataframe with columns Type, Tenor and Rate with Type in (LIBOR, SWAP), Tenor the instrument maturity in years and rate the rate per annum
#' 
fGetCashflowsLibor <- function( dfInstrument ) {
	
	dcf <- dfInstrument[["Tenor"]] 
	cashflows <- 1 + dfInstrument[["Rate"]] * dcf
	times <- fGetTimesLibor( dfInstrument )
	return( data.frame( times, cashflows ) )
	
}


#' Gets the cashflow schedule for a swap
#' 
#' @param dfInstrument A set of market instruments as a dataframe with columns Type, Tenor and Rate with Type in (LIBOR, SWAP), Tenor the instrument maturity in years and rate the rate per annum
#' 
fGetCashflowsSwap <- function( dfInstrument ) {
	
	freq <- dfInstrument[["Frequency"]]
	dcf <- 1 / dfInstrument[["Frequency"]]
	cashflows <- rep( dfInstrument[["Rate"]] * dcf, dfInstrument[["Tenor"]] * freq )
	cashflows[ length( cashflows ) ] <- cashflows[ length( cashflows ) ] + 1
	times <- fGetTimesSwap( dfInstrument )
	
	if ( ( 0L == length( cashflows ) ) || ( 0L == length( times ) ) ) stop( "No cashflows calculated for swap, check Tenor and Frequency" )
	
	return( data.frame( times, cashflows ) )
	
}
