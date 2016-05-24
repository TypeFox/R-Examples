# TODO: Add comment
# 
# Author: ecor
###############################################################################


NULL

#'
#' This function saves all spatially distrubuted information contained in the recovery folder into a comprehensive \code{list} object.
#' 
#' 
#' 
#' 
#' 
#' 
#' @author Emanuele Cordano
#' 
#' @param recFolder directory when recvery maps are set. In GEOtop it is ...
#' @param xx charcter String. Default is \code{"0000"}
#' @param extension file estension used for asccii recovery map files. It must contains \code{'.'} as the first character. Defaut is \code{".asc"} . 
#' @param formatter string character for the the decimal formatter to be used. Default is \code{"L\%04d"}.
#' @param nsoillayers number of soil layers used in the GEOtop simulation
#' @param ... further arguments 
#' 
#' @return a \code{list} object containining all recovery raster maps. 
#' 
#' @note This function has been used with the built 1.225-9 of GEOtop . 
#' @seealso \code{\link{brick.decimal.formatter}},
#' 
#' \code{\link{raster}},\code{\link{set.geotop.recovery.state}},
#' 
#' \code{\link{write.vectorized.geotop.recovery}},\code{\link{read.vectorized.geotop.recovery}}
#' 
#' @export
#' @examples
#' library(geotopbricks)
#' example_Rscript <- system.file('template/example.geotop.recovery.state.R',package="geotopbricks")
#' example_Rscript
#' 
#' # Not Run because it elapses too long time!!! 
#' # Please Uncomment the following line to run by yourself!!!
#' # source(example_Rscript)
#' 

 




get.geotop.recovery.state <- function(recFolder,xx="0000",formatter="L%04d",extension=".asc",nsoillayers=10,...) {


	names<-c("RainOnCanopy","SnowAge","SnowIceContent","SnowLayersNumber","SnowLiqWaterContent","SnowOnCanopy","SnowTemperature","SnowThickness","SoilChannelIceContent","SoilChannelPressure","SoilChannelTemperature","SoilIceContent","SoilPressure","SoilTemperature","VegTemperature")      
#	xx <- "0000"
    layer <- array(formatter,length(names))
	
	# 4 Groups of rasterbricks or rasters !!! 
	noLayers <- (names %in% c("RainOnCanopy","VegTemperature","SnowAge","SnowLayersNumber","SnowOnCanopy"))	
	soilLayersWithZero <- (names %in% c("SoilPressure","SoilChannelPressure"))
	soilLayers <- (str_detect(names,"Soil")) & (!soilLayersWithZero) & (!noLayers)
	snowLayers <- str_detect(names,"Snow") & (!noLayers) & (!soilLayers) & (!soilLayersWithZero)
	
	layer[noLayers] <- ""
	files <- paste(names,xx,layer,extension,sep="")
	files_w <- paste(recFolder,files,sep="/")
#	## ##print(files)
	out <- list()
	
	for (it in names[noLayers]) {
		
		

		x <- as.character(files_w[names==it])
		 print(x)
		out[it] <- raster(x)
		
		
	}
	for (it in names[soilLayersWithZero]) {
		
		x <- as.character(files_w[names==it])
		
		out[it] <- brick.decimal.formatter(file=x,nlayers=nsoillayers,start.from.zero=TRUE)
	
	}
	for (it in names[soilLayers]) {
	
		x <- as.character(files_w[names==it])
		
		out[it] <- brick.decimal.formatter(file=x,nlayers=nsoillayers,start.from.zero=FALSE)
	
	}
	
	# number of snow leyers is detected by raster layer 'out$SnowLayersNumber'
	out$SnowLayersNumber <- setMinMax(out$SnowLayersNumber)
	nsnowlayers <- maxValue(out$SnowLayersNumber)
	nsnowlayers[nsnowlayers<1] <- 1
	
	for (it in names[snowLayers]) {
		

		x <- as.character(files_w[names==it])
	
		out[it] <- brick.decimal.formatter(file=x,nlayers=nsnowlayers,start.from.zero=FALSE)
		
		
		
	}
	
	out$names <- names
	out$files <- files
	
	out$noLayers <- noLayers
	out$soilLayersWithZero <- soilLayersWithZero
	out$soilLayers <- soilLayers
	out$snowLayers <- snowLayers
	
	
	
	return(out)
	
}
