#' Aligns raster files
#' 
#' Aligns a raster to a reference raster.
#' 
#' Aligns a raster to the extent and projection of a reference raster and 
#' matches the resolution of the reference raster.  This is helpful in 
#' preparing multiple files of different projections, resolutions, extents,
#' and rotations for performing map algebra or change detection.
#' 
#' @param unaligned Character.  The filename of a raster to be aligned to the reference raster.
#' @param reference Character. The filename of a raster to be used as the reference for the alignment.
#' Syncing will use the reference's projection, resolution, and extent.
#' @param dstfile Character. The filename of the synchronized output file.
#' @param output_Raster Logical. Return output dst_dataset as a RasterBrick?
#' @param nThreads Numeric or Character. If numeric, the number of threads to use.  Setting to > 1 enables multithreaded execution.  Can also be "ALL_CPUS" to use all available CPUS. Default is 1.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @param ... parameters to be passed to {\link{gdalwarp}} (e.g. resampling approach).
#' @return Either NULL or a RasterBricks depending on whether output_Raster is set to TRUE.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net})
#' @seealso \code{\link{gdalwarp}}
#' @import raster
#' @export

align_rasters <- function(unaligned,reference,dstfile,output_Raster=FALSE,
		nThreads=1,
		verbose=FALSE,...)
{
	# Get projection from reference
	reference_info <- gdalinfo(reference,proj4=TRUE,raw_output=FALSE,verbose=verbose)
	
#	cmd_output <- gdalinfo(reference,verbose=TRUE,proj4=TRUE,raw_output=FALSE)
	
	proj4_string <- reference_info$proj4
	
	bbox <- reference_info$bbox
	te <- c(reference_info$bbox[1,1],reference_info$bbox[2,1],reference_info$bbox[1,2],reference_info$bbox[2,2])
	ts <- c(reference_info$columns,reference_info$rows)
	
	if(missing(dstfile))
		dstfile <- tempfile()
	
	if(is.character(nThreads))
	{
		if(nThreads=="ALL_CPUS")
		{
			multi=TRUE
			wo="NUM_THREADS=ALL_CPUS"
		}
	} else
	{
		if(nThreads==1)
		{
			multi=FALSE
			wo=NULL
		}
		else
		{
			multi=TRUE
			wo=paste("NUM_THREADS=",nThreads,sep="")
		}
	}
	
	synced <- gdalwarp(srcfile=unaligned,dstfile=dstfile,te=te,t_srs=proj4_string,
			ts=ts,output_Raster=output_Raster,multi=multi,wo=wo,verbose=verbose,...)
	return(synced)
	
}