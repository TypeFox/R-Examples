#' Construct the Smith-Wilson yield curve
#' 
#' Constructs the SmithWilson ZCB function based on the given market inputs and parameter choices
#' 
#' @return a list containing:
##' \itemize{
##'  \item{"P"}{ a function of time which gives the ZCB price to that term }
##'  \item{"xi"}{ the vector of weights applied to the kernel functions to obtain the ZCB price }
##'  \item{"K"}{ the (compound) kernel vector }
##' } 
#' 
#' @param TimesVector A vector of all cashflow times
#' @param CashflowMatrix A matrix of all cashflows, instruments in rows, times in columns
#' @param MarketValueVector A vector of market values of the insturments
#' @param ufr The Ultimate Forward Rate (UFR) of the Smith-Wilson kernel 
#' @param alpha The rate of reversion of forward rates to the UFR in the Smith-Wilson kernel
#' 
#' @export
#' 
fFitSmithWilsonYieldCurve <- function( TimesVector, CashflowMatrix, MarketValueVector, ufr, alpha ) {
	
	if ( ufr < 0 ) warning( "Parameter ufr should be greater than zero" )
	if ( alpha <= 0 ) stop( "Parameter alpha must be greater than zero" )
	
	fBase <- function( t ) { exp( -ufr * t) }
	fKernel <- function(t, u) { fWilson( t, u, ufr, alpha ) }
		
	SmithWilsonYieldCurve <- fFitYieldCurve( TimesVector, CashflowMatrix, MarketValueVector, fKernel, fBase )
	
	class( SmithWilsonYieldCurve ) <- c( "SmithWilsonYieldCurve", "YieldCurve" )
	
	return( SmithWilsonYieldCurve )
	
}