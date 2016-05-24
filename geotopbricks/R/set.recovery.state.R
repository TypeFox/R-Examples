# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL

#'
#' This function re-writes the recovery ascii raster maps in a given folder 
#' 
#' 
#' 
#' 
#' 
#' 
#' @author Emanuele Cordano
#' 
#' @param rec a \code{list} object returened by \code{\link{get.geotop.recovery.state}}
#' @param newRecFolder directory where to write all recovery raster asccii maps
#' @param ... further arguments 
#' 
#' @seealso \code{\link{get.geotop.recovery.state}},\code{\link{writeRasterxGEOtop}}
#' @export
#' @examples
#' # See the examples of the 'get.geotop.recovery.state' function


set.geotop.recovery.state <- function(rec,newRecFolder,...) {
		out <- 0 
		NoLayers <- rec$noLayers 
		LayersWithoutZero <-  rec$soilLayers | rec$snowLayers
		LayersWithZero <- rec$soilLayersWithZero
		
		names <- rec$names
		files_w <- paste(newRecFolder,rec$files,sep="/")

		for (it in names[NoLayers]) {
			
			file  <- files_w[names==it]

			####print(it)
			writeRasterxGEOtop(x=rec[[it]],filename=file,use.decimal.formatter=FALSE,start.from.zero=FALSE,...)								
			
		}
		
		
		for (it in names[LayersWithoutZero]) {
			
			file  <- files_w[names==it]
			## print(it)
			## print(file)
			writeRasterxGEOtop(x=rec[[it]],filename=file,use.decimal.formatter=TRUE,start.from.zero=FALSE,...)								
			
		}
		
		for (it in names[LayersWithZero]) {
			
			file  <- files_w[names==it]	
			### print(it)
			## print(file)
			
			writeRasterxGEOtop(x=rec[[it]],filename=file,use.decimal.formatter=TRUE,start.from.zero=TRUE,...)
			
		}		
#	out$names <- names
#	out$files <- files
	
#	out$noLayers <- noLayers
#	out$soilLayersWithZero <- soilLayersWithZero
#	out$soilLayers <- soilLayers
#	out$snowLayers <- snowLayers
	
	
	
	return(out)
	
}
