#' Constructs the ZCB function based on the given market inputs and a specific kernel and base function
#' 
#' 
#' @param TimesVector A vector of all cashflow times
#' @param CashflowMatrix A matrix of all cashflows, instruments in rows, times in columns
#' @param MarketValueVector A vector of market values of the insturments
#' @param fKernel a function of two times used as the Kernel "basis" function
#' @param fBase a function giving the base level of the curve
#' 
#' @return a list comprising elements: a function of time which gives the ZCB price to that time
#' 
fFitYieldCurve <- function( TimesVector, CashflowMatrix, MarketValueVector, fKernel, fBase ) {
	
	KernelFunctionMatrix <- fCreateKernelMatrix( TimesVector, fKernel )
	
	BaseZeroVector <- fBase( TimesVector )
	KernelWeights <- fFitKernelWeights( CashflowMatrix, KernelFunctionMatrix, MarketValueVector, BaseZeroVector)
	
	fCompoundKernel <- function( t ) { CashflowMatrix %*% fKernel( t, TimesVector ) }
	
	fPricingFunction <- function( t ) { fBase( t ) + t( KernelWeights ) %*% fCompoundKernel( t ) }
	
	YieldCurve <- list()
	YieldCurve$P <- fPricingFunction
	YieldCurve$xi <- KernelWeights
	YieldCurve$K <- fCompoundKernel

	return( YieldCurve )
}