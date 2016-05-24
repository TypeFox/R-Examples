NULL 

#'
#' Graphic Representation of a Color legend of a Raster or GeotopbrickRaster object as a Color bar, inspired by the function written by John Colby
#' 
#' @param x a Raster or GeotopRasterBrick object
#' @param col the color palette used 
#' @param max,min maximum and minimum value (used if you need to crop the scale legend within a cartein interval)
#' @param ...  arguments to be passed to \code{\link{color.bar}}
#' 
#' @seealso \code{\link{color.bar}},\code{\link{setMinMax}}
#' @export 
#' 
#' @examples
#' library(geotopbricks)
#' 
#' ## Simulation working path
#'
#' file <- system.file("rendena100/SnowDepthMapFile-2014-MA-mean-winter-2013-2014.asc",
#' package="geotopbricks")
#' snow <- raster(file)
#' 
#' min <- 0 # snow depth expressed in millimeters
#' max <- 2500 # snow depth expressed in millimeters
#' 
#' colors <- terrain.colors(1000)
#' 
#' color.bar.raster(x=snow,col=colors,digits=2)
#' color.bar.raster(x=snow,col=colors,min=min,max=max,digits=2)
#' 


color.bar.raster <- function(x,col,min=NA,max=NA,...) {
	
	smin <- min
	smax <- max
	
	if (class(x)=="GeotopRasterBrick") {
	
		min <- min_value(x)
		max <- max_value(x)
	
	} else {
		x <- setMinMax(x)		
		min <- min(minValue(x),na.rm=TRUE)	
		max <- max(maxValue(x),na.rm=TRUE)
		
	}
	
	if (!is.na(smin)) {
		
		v <- seq(from=min,to=max,length.out=length(col))
		
		index <- which(v>=smin)
		
		col <- col[index]
		v <- v[index]
		min <- min(v,na.rm=TRUE)
	
	}
	if (!is.na(smax)) {
		
		v <- seq(from=min,to=max,length.out=length(col))
		
		index <- which(v<=smax)
		
		col <- col[index]
		v <- v[index]
		max <- max(v,na.rm=TRUE)
		
	}
	
	
	color.bar(lut=col,min=min,max=max,...)
	
	# TO DO 
}