NULL


#'  
#' Extracts generated time series of Daily Minimum Temperature from a random multi-realization obtained by \code{\link{generateTemperatureTimeseries}} function 
#' 
#' @param res_multigen matrix containing standardized values of daily temperature as returned by \code{\link{generateTemperatureTimeseries}} (first item)
#' @param std vector containing standard deviation for each minimun temperature anomalies
#' @param SplineAdv matrix containing the averaged daily values of minimum temperature obtained by a spline interpolation of the monthly climate
#' 
#' 
#' @export 
#' 
#' @author  Emanuele Cordano, Emanuele Eccel
#'    
#'   
#'
#'        
#' 
#' @return  a matrix with generated minimum temperature




extractTnFromAnomalies <-
function (res_multigen,std,SplineAdv) {

	
	ntall <- as.integer(ncol(res_multigen))
	ntn <- as.integer(ncol(res_multigen)/2)
	
	nstd <- length(std)
	out <- res_multigen[,(ntn+1):ntall]%*% diag(x=std,nrow=nstd) + SplineAdv
	return(out)
	
}

