#' Extract filenames from all Raster* objects.
#' 
#' @param x Raster*.  A Raster* object (even one without values/in memory) to determine the filename(s).
#' @param unique Logical. Only return unique filenames?  If FALSE, one filename per layer.
#' 
#' @return Character vector of filenames.
#' @author Jonathan A. Greenberg
#' @seealso \code{\link[raster]{filename}}
#' @details This is an expansion of filename() that allows for RasterStacks, in-memory Raster*s,
#' and Raster*s without values.  If a filename is not found, the entry will be "".
#' 
#' @examples { 
#' tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' raster_to_filenames(tahoe_highrez)
#' raster_to_filenames(tahoe_highrez,unique=TRUE)
#' nodata <- raster()
#' raster_to_filenames(nodata)
#' }
#' @import raster
#' @export

raster_to_filenames <- function(x,unique=FALSE)
{
	if(!hasValues(x)) return("")
	if(inMemory(x)) return("")
	
	filenames <- sapply(X=seq(nlayers(x)),
			function(X,raster)
			{
				raster_layer <- raster(raster,layer=X)
				if(!hasValues(raster_layer)) return("")
				if(inMemory(raster_layer)) return("")
				else
				return(filename(raster_layer))
			},			
			raster=x)
	
	if(unique) 
	{
		filenames <- unique(filenames)
	} 
	return(filenames)
}