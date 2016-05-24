#' Mosaic raster files using GDAL Utilities
#' 
#' @param gdalfile Character. Input files (as a character vector) or a wildcard search term (e.g. "*.tif") 
#' @param dst_dataset Character. The destination file name.
#' @param output.vrt Character. Output VRT file.  If NULL a temporary .vrt file will be created.
#' @param output_Raster Logical. Return output dst_dataset as a RasterBrick?
#' @param trim_margins Numeric. Pre-crop the input tiles by a fixed number of pixels before mosaicking.  Can be a single value or four values representing the left, top, right, and bottom margins, respectively.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @param ... Parameters to pass to \code{\link{gdalbuildvrt}} or \code{\link{gdal_translate}}.
#' 
#' @details This function mosaics a set of input rasters (gdalfile) using parameters
#' found in \code{\link{gdalbuildvrt}} and subsequently exports the mosaic to 
#' an output file (dst_dataset) using parameters found in \code{\link{gdal_translate}}.  The user
#' can choose to preserve the intermediate output.vrt file, but in general this is not
#' needed.
#' @return Either a list of NULLs or a list of RasterBricks depending on whether output_Raster is set to TRUE.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net})
#' @seealso \code{\link{gdalbuildvrt}}, \code{\link{gdal_translate}}
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(raster) && require(rgdal) && valid_install)
#' {
#' layer1 <- system.file("external/tahoe_lidar_bareearth.tif", package="gdalUtils")
#' layer2 <- system.file("external/tahoe_lidar_highesthit.tif", package="gdalUtils")
#' mosaic_rasters(gdalfile=c(layer1,layer2),dst_dataset="test_mosaic.envi",separate=TRUE,of="ENVI",
#' 		verbose=TRUE)
#' gdalinfo("test_mosaic.envi")
#' }
#' @import rgdal
#' @importFrom utils write.table
#' @export

mosaic_rasters <- function(gdalfile,dst_dataset,output.vrt=NULL,output_Raster=FALSE,
		trim_margins = NULL,
		verbose=FALSE,
		...)
{
	# CRAN check to fix foreach variable errors:
	k <- NULL	
	
	# Check to make sure all the input files exist on the system:
	if(verbose) message("Checking to make sure all the input files exist...")
	files_exist <- sapply(gdalfile,file.exists)
	if(!all(files_exist))
	{
		missing_files <- gdalfile[!files_exist]
		stop(paste("Some of the input files are missing:",missing_files))
	}
	
	if(output_Raster && (!requireNamespace("raster") || !requireNamespace("rgdal")))
	{
		warning("rgdal and/or raster not installed. Please install.packages(c('rgdal','raster')) or set output_Raster=FALSE")
		return(NULL)
	}
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation()
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	if(is.null(output.vrt))
	{
		output.vrt <- paste(tempfile(),".vrt",sep="")
	}

	# Shrink 
	if(!is.null(trim_margins))
	{
		if(verbose) message("Trimming margins...")
		if(length(trim_margins)==1) 
		{
			trim_margins <- rep(trim_margins,4)
		}
		
		gdalfile_vrts <- foreach(k = seq(gdalfile),.combine="c") %dopar%
				{
					temp_gdalfile_info <- gdalinfo(gdalfile[k],raw_output=F)
					new_xmin <- temp_gdalfile_info$bbox[1,1] + trim_margins[1]*abs(temp_gdalfile_info$res.x)
					new_xmax <- temp_gdalfile_info$bbox[1,2] - trim_margins[3]*abs(temp_gdalfile_info$res.x)
					new_ymin <- temp_gdalfile_info$bbox[2,1] + trim_margins[2]*abs(temp_gdalfile_info$res.y)
					new_ymax <- temp_gdalfile_info$bbox[2,2] - trim_margins[4]*abs(temp_gdalfile_info$res.y)
					
					temp_vrt_name <- paste(tempfile(),".vrt",sep="")
					
					gdalbuildvrt(gdalfile=gdalfile[k],
							output.vrt=temp_vrt_name,te=c(new_xmin,new_ymin,new_xmax,new_ymax),
							verbose=verbose)
					
					return(temp_vrt_name)
				}
		gdalfile <- gdalfile_vrts
		
	}
	
	# There is an error that occurs with a lot of files.  We are going to fix this by
	#	creating a file list.
	
	temp_file_list_name <- paste(tempfile(),".txt",sep="")
	write.table(gdalfile,temp_file_list_name,row.names=F,col.names=F,quote=F)
	
	# Now pass the right arguments to each function:
#browser()
#	additional_arguments <- list(...)
#	gdalbuildvrt_formals <- names(formals(gdalbuildvrt))
#	gdal_translate_formals <- names(formals(gdal_translate))
#	
#	gdalbuildvrt_additional_args <- additional_arguments[names(additional_arguments) %in% gdalbuildvrt_formals]
#	gdal_translate_additional_args <- additional_arguments[names(additional_arguments) %in% gdal_translate_formals]
#	
#	
	gdalbuildvrt(input_file_list=temp_file_list_name,output.vrt=output.vrt,verbose=verbose,...)
	outmosaic <- gdal_translate(src_dataset=output.vrt,dst_dataset=dst_dataset,
		output_Raster=output_Raster,verbose=verbose,...)
	return(outmosaic)
	
}