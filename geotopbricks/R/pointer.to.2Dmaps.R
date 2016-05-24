NULL
#'
#'@param wpath complete working path to *.asc maps are saved 
#'@param map.prefix string prefix name map before 
#'@param suffix z-time or time suffix plus file extention character string. Default for GEOtop application is \code{"L\%04dN\%04d.asc"} for xy+z+time maps or \code{"N\%04d.asc"}  for xy+time maps.
#'@param zoo.index time or date index. Default is \code{NULL} , otherwise function returns a zoo object with \code{zoo.index} as index.
#'@param ntime number of time instant. If \code{zoo.index} is not \code{NULL}, it is calculated from \code{zoo.index} length.
#'@param nlayers number of vertical layers. 
#'
#'@author Emanuele Cordano
#'
#' @return A dota.frame or zoo object containig the paths to maps fpr each time and z layer.
#' @title pointer.to.maps.xyz.time
#' @name pointer.to.maps.xyz.time
#' @rdname pointer.to.maps.xyz.time
#' @export 
#'

pointer.to.maps.xyz.time <- function(wpath,map.prefix="thetaliq",suffix="L%04dN%04d.asc",zoo.index=NULL,ntime,nlayers) {
	
	if (!is.null(zoo.index)) ntime <- length(zoo.index)
	
	out <- (as.data.frame(
						array(paste(wpath,"/",map.prefix,suffix,sep=""),c(ntime,nlayers))
				))
	

	
	for (l in 1:ncol(out)) {
		
		out[,l] <- as.character(out[,l])
		
	}
	
	for (l in 1:ncol(out)) {
		for (n in 1:nrow(out)) {
			
			out[n,l] <- sprintf(as.character(out[n,l]),l,n)	
		}
		
	} 
	
	if (!is.null(zoo.index)) {
		
		out <- as.zoo(out)
		index(out) <- zoo.index
	}
	
	
	
	
	return(out)
	
}

#'@title pointer.to.maps.xy.time
#' @name pointer.to.maps.xy.time
#' @rdname pointer.to.maps.xyz.time
#' @export
NULL 



pointer.to.maps.xy.time <- function(wpath,map.prefix="SWE",suffix="N%04d.asc",zoo.index=NULL,ntime) {
	
	if (!is.null(zoo.index)) ntime <- length(zoo.index)
	
	out <- array(paste(wpath,"/",map.prefix,suffix,sep=""),c(ntime))			
	
	for (n in 1:length(out)) {
		
		out[n] <- as.character(out[n])
		out[n] <- sprintf(out[n],n)
	}
	
	
	if (!is.null(zoo.index)) {
		
		out <- as.zoo(out)
		index(out) <- zoo.index
	}
	
	
	
	
	return(out)
	
}







