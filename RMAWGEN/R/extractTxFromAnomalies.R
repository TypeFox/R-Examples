NULL


#'  
#' Extracts generated time series of Daily Maximum Temperature from a random multi-realization obtained by \code{\link{generateTemperatureTimeseries}} function 
#' 
#' @param res_multigen matrix containing standardized values of daily temperature as returned by \code{\link{generateTemperatureTimeseries}} (first item)
#' @param std vector containing standard deviation for each maximum temperature anomalies
#' @param SplineAdv matrix containing the averaged values of maximum temperature obtained by a spline interpolation of monthly climate
#' 
#' @export 
#' 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'    
#'   
#'
#'
#'        
#' @return a matrix with generated maximum temperature





extractTxFromAnomalies <-
function (res_multigen,std,SplineAdv) {

	
	ntall=as.integer(ncol(res_multigen))
	ntn=as.integer(ncol(res_multigen)/2)
	
	nstd <- length(std)

	out <- res_multigen[,1:ntn]%*% diag(x=std,nrow=nstd) + SplineAdv
	return(out)
	
}

