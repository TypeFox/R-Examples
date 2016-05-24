#' Create the matrix of kernel functions
#' 
#' Creates a J x J matrix [ w(u_i,u_j) ] where J is the number of cashflow times in the calibration set
#' 
#' @param times a vector of cashflow times
#' @param fKernel a kernel to apply (a function of times x times returning a matrix )
#' 
fCreateKernelMatrix <- function( times, fKernel ) {
	
	W <- fKernel( times, times )
	
	return( W )
	
}