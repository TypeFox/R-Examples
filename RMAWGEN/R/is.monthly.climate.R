
NULL


 
#' 
#' Verifies if 'climate' represents the monthly climatology in one year, i.e  'climate' is monthly.climate type matrix whose rows represent months and each column represents a station. It is also used in \code{\link{setComprehensiveTemperatureGeneratorParameters}}.
#' 
#'  
#' @param  climate matrix containing the 'monthly climatology' data 
#' @param  nstation number of variable measurement stations (columns of the matrix 'climate')
#' @param  nmonth number of months in one year (it can be different if climate is represented by seasonal avarages or others), Default is 12 (recommended). (it can be different if climate is represented by seasonal averages, in this case 4)  
#' @param  verbose Prints output and warining messagrs only if is \code{TRUE}.
#'    
#' @export 
#'       
#' @author  Emanuele Cordano, Emanuele Eccel  
#' 
#' @return  A logical variable if the matrix 'climate' is monthly.climate type
#' 
#' @seealso \code{\link{setComprehensiveTemperatureGeneratorParameters}}
# @examples 
# data(trentino_predictions)
# is.monthly.climate(Tn_2021_2050_50,nstation=ncol(Tn_2021_2050_50))
#



is.monthly.climate <-
function (climate,nstation=3,nmonth=12,verbose=TRUE) {

	
	
	if (is.null(climate)) {
		return(FALSE)
	} else if (!is.matrix(climate)) {
		if (is.list(climate)) {
			
			vec <- array(FALSE,length(climate))
			if (is.matrix(vec[[1]])) {
				for (i in 1:length(vec)) {
					vec[i] <- is.monthly.climate(vec[[i]],nstation=nstation,nmonth=nmonth,verbose=verbose)
				}
				return(as.logical(max(vec)))			
			} else {
				if (verbose) print("Error: The format for monthly mean climate is not a matrix !!")
			}
			
			
			
		} else {
			if (verbose) print("Error: The format for monthly mean climate is not a matrix !!")
			return(FALSE)
		}
	} else {
		v <- dim(climate)
		if ((v[1]!=nmonth) | (v[2]!=nstation)) {
			if (verbose) print("Error: The matrix format for monthly mean climate is not correct !!")
			return(FALSE)
		} else {
			return(TRUE)
		}
		
		return(FALSE)
	}
	
	
	return(FALSE)
	
}

