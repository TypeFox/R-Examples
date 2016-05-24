#' Performs an area-weighted resampling of raster datasets.
#' 
#' @param from Raster* The sources raster to be resampled.
#' @param to Raster* A target raster that the from will be resampled to (extent, resolution, projection).
#' @param method Character. Default is "mode". See details.
#' @param na.rm Logical. Remove NAs before calculating cell stats?
#' @param verbose logical. Enable verbose execution? Default is FALSE.  
#' @param ... Currently unsupported.
#' @author Jonathan A. Greenberg
#' @seealso \code{\link[raster]{projectRaster}}, \code{\link[raster]{extract}}, \code{\link[raster]{aggregate}}
#' @details This function is designed to solve the problem of resampling/reprojecting rasters
#' using area-based, not point based (e.g. nearest neighbor, bilinear, cubic convolution), 
#' resampling.  The output pixel is a function of the areas of the input pixels, so this
#' should be used for resampling from a finer resolution to a coarser resolution.
#' 
#' The method defaults to "mode", which will return the value covering the largest area
#' of the output pixel area.  Other methods will be added in the future.
#' 
#' A word of warning: this algorithm is SLOW.  The function uses focal_hpc, so we 
#' highly recommend using it with a foreach engine running (e.g. use sfQuickInit() ).
#' Keep in mind this is a "dirty" parallel problem, so different chunks may execute at
#' different speeds and have different memory footprints.
#' 
#' @import raster
#' @import foreach
#' @export

projectRaster_rigorous <- function(from,to,method="mode",na.rm=FALSE,verbose=FALSE,...)
{
	chunk_function <- function(x,from,method,na.rm,...)
	{
#		gc()
		if(class(from)!="RasterLayer") stop("Currently, from= must be a RasterLayer.")
		
		# This should be changed and the polys made "manually":
		chunk_vector <- rasterToPolygons(x,na.rm=FALSE,n=16)
		chunk_vector_reproject <- spTransform(chunk_vector,CRS(projection(from)))
		chunk_vector_extract <- extract(from,chunk_vector_reproject,weights=TRUE,na.rm=FALSE)
		if(method=="mode")
		{
			chunk_vector_extract_area <- 
					sapply(chunk_vector_extract,
							function(x)
							{
								sum_class<-tapply(x[,2],x[,1],sum)
								class_names <- names(sum_class)
								return(as.numeric(class_names[which.max(sum_class)]))
							}
					)
		}
		
		if(method=="mean")
		{
			chunk_vector_extract_area <- 
					sapply(chunk_vector_extract,
							function(x,na.rm)
							{
								return(weighted.mean(x[,1],x[,2],na.rm))
#								sum_class<-tapply(x[,2],x[,1],sum)
#								class_names <- names(sum_class)
#								return(as.numeric(class_names[which.max(sum_class)]))
							},na.rm=na.rm
					)
		}
		
		return(array(chunk_vector_extract_area,dim=c(dim(x)[2],dim(x)[1],1)))
	}
	if(class(to)!="RasterLayer")
	{
		x <- raster(to,layer=1)
	} else
	{
		x <- to
	}
	
	return(focal_hpc(x=x,fun=chunk_function,args=list(from=from,method=method,na.rm=na.rm),
					chunk_format="raster",blocksize=1,verbose=verbose,...))
}
