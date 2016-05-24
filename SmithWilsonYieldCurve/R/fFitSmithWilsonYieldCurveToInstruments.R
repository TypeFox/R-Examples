#' Construct the Smith-Wilson yield curve
#' 
#' Constructs the SmithWilson ZCB function based on the given market inputs and parameter choices.
#' Primarily a convenience wrapper around other package functions
#' 
#' @return a list containing:
##' \itemize{
##'  \item{"P"}{ a function of time which gives the ZCB price to that term }
##'  \item{"xi"}{ the vector of weights applied to the kernel functions to obtain the ZCB price }
##'  \item{"K"}{ the (compound) kernel vector }
##' } 
#' 
#' @param InstrumentSet A set of market instruments as a dataframe with columns
#' 		 \itemize{
#'	 		 \item{"Type"}{One of (LIBOR, SWAP) }
#' 			 \item{"Tenor"}{The instrument maturity in years}
#' 			 \item{"Frequency"}{The payment frequency (ignored for Type=="LIBOR" )}
#' 		 	\item{"Rate"}{The coupon rate per annum in percent}
#' 		 }
#' @param ufr The Ultimate Forward Rate (UFR) of the Smith-Wilson kernel 
#' @param alpha The rate of reversion of forward rates to the UFR in the Smith-Wilson kernel
#' 
#' @export
#' 
fFitSmithWilsonYieldCurveToInstruments <- function( InstrumentSet, ufr, alpha ) {

	MarketValueVector <- rep(1, length.out=NROW(InstrumentSet) )
	CashflowMatrix <- fCreateCashflowMatrix( InstrumentSet )
	TimesVector <- fCreateTimeVector( InstrumentSet )
	
	SmithWilsonYieldCurve <- fFitSmithWilsonYieldCurve( TimesVector, CashflowMatrix, MarketValueVector, ufr, alpha )
	
	return( SmithWilsonYieldCurve )
	
}