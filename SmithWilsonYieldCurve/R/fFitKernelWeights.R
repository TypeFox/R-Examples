#' Solve for the vector xi of kernel weights
#' 
#' @param CashflowMatrix A matrix of all cashflows, instruments in rows, times in columns
#' @param KernelFunctionMatrix A matrix of kernel function values
#' @param MarketValueVector A vector of market values of the insturments
#' @param BaseZeroVector A vector of "base" values for the zeros
#' 
fFitKernelWeights <- function( CashflowMatrix, KernelFunctionMatrix, MarketValueVector, BaseZeroVector ) {
	
	xi <- solve( CashflowMatrix %*% KernelFunctionMatrix %*% t( CashflowMatrix ) ) %*% ( MarketValueVector - CashflowMatrix %*% BaseZeroVector )
	
	return( xi )
	
}