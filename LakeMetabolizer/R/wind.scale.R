#'@title Wind Scaling U10 - exponential conversion to 10m wind speed
#'@description Scale wind speed to standard U10 (10 meters) based on height of observations
#'@name wind.scale
#'@aliases
#'wind.scale
#'wind.scale.base
#'@usage
#'## Used for timeseries data in a data.frame 
#'wind.scale(ts.data, wnd.z)
#'
#'## Used for raw numeric data
#'wind.scale.base(wnd, wnd.z)
#'@param ts.data Object of class data.frame containing a wnd column.
#'@param wnd measured wind speed (Units: typically m s-1, but it is unit agnostic)
#'@param wnd.z height of anemometer (Units: meters)
#'@return
#'## wind.scale
#'Returns a data frame with columns datetime and wnd_10 and the same number of rows as ts.data
#'
#'## wind.scale.base
#'Returns a vector with the same length as wnd
#'@author
#'Aline Jaimes, Luke A. Winslow
#'@references
#'Saucier, W. 2003. \emph{Principles of Meteorological Analysis}. Dover Publications. New York. p433
#'@seealso
#'Models of gas flux \link{k.cole}, \link{k.crusius}, \link{k.macIntyre}, & \link{k.read}.
#'@details This function transforms wind speed to the standard U10, 
#'speed at 10 meters, based on the common exponential wind profile assumption. 
#'wind.scale defaults to using the supplied wnd.z value. If wnd.z is not supplied, 
#'it attempts to determine the anemometer height from the suffix of the header 
#'(e.g., a header of wnd_3 would mean an anemometer height of 3 meters).
#'@examples
#'wndSpeed <- c(5.1,6.3,6.3,5.2,7,7.2)
#'wndHeight <- 2
#'
#'wind.scale.base(wndSpeed, wndHeight)
#'@export
wind.scale <- function(ts.data, wnd.z){
	wnd = get.vars(ts.data, 'wnd')
	
	if(missing(wnd.z)){
		wnd.z = get.offsets(wnd)
    if(is.na(wnd.z)){
			stop('Unknown wind height. Must supply wnd.z parameter or have offset defined in header of ts.data.')
		}
	}
	
	if(ncol(wnd) > 2 || length(wnd.z) > 1){
		stop('too many columns supplied to scale.exp.wind. Please supply only one datetime and wnd columns.')
	}
	
	u10 = wind.scale.base(wnd[,2], wnd.z)
	
	return(data.frame(datetime=ts.data$datetime, wnd_10=u10))
}

#Reference for this
#Arya 1988 (Introduction to micrometeorology)
#'@export
wind.scale.base <- function(wnd, wnd.z){
	U10 <- wnd * (10/wnd.z)^(0.15) # I'm pretty sure that (1/7) should be 0.15.  (1/7) used to be 1.7, but I think it was incorrectly copied from the U10^1.7 of k.cole.R. I have other code that has this value as 0.15.
	return(U10)
}