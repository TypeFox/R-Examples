#' gdalbuildvrt
#' 
#' R wrapper for gdalbuildvrt: Builds a VRT from a list of datasets
#' 
#' @param gdalfile Character. Input files (as a character vector) or a wildcard search term (e.g. "*.tif") 
#' @param output.vrt Character. Output VRT file.
#' @param tileindex Logical. Use the specified value as the tile index field, instead of the default value with is 'location'.
#' @param resolution Character. ("highest"|"lowest"|"average"|"user") In case the resolution of all input files is not the same, the -resolution flag enables the user to control the way the output resolution is computed. 'average' is the default. 'highest' will pick the smallest values of pixel dimensions within the set of source rasters. 'lowest' will pick the largest values of pixel dimensions within the set of source rasters. 'average' will compute an average of pixel dimensions within the set of source rasters. 'user' is new in GDAL 1.7.0 and must be used in combination with the -tr option to specify the target resolution.
#' @param te Numeric. c(xmin,ymin,xmax,ymax) (starting with GDAL 1.7.0) set georeferenced extents of VRT file. The values must be expressed in georeferenced units. If not specified, the extent of the VRT is the minimum bounding box of the set of source rasters.
#' @param tr Numeric. c(xres,yres) (starting with GDAL 1.7.0) set target resolution. The values must be expressed in georeferenced units. Both must be positive values. Specifying those values is of course incompatible with highest|lowest|average values for -resolution option.
#' @param tap Logical. (GDAL >= 1.8.0) (target aligned pixels) align the coordinates of the extent of the output file to the values of the -tr, such that the aligned extent includes the minimum extent.
#' @param separate Logical. (starting with GDAL 1.7.0) Place each input file into a separate stacked band. In that case, only the first band of each dataset will be placed into a new band. Contrary to the default mode, it is not required that all bands have the same datatype.
#' @param b Numeric. (GDAL >= 1.10.0) Select an input band to be processed. Bands are numbered from 1. If input bands not set all bands will be added to vrt
#' @param sd Numeric. (GDAL >= 1.10.0) If the input dataset contains several subdatasets use a subdataset with the specified number (starting from 1). This is an alternative of giving the full subdataset name as an input.
#' @param allow_projection_difference Logical. (starting with GDAL 1.7.0) When this option is specified, the utility will accept to make a VRT even if the input datasets have not the same projection. Note: this does not mean that they will be reprojected. Their projection will just be ignored.
#' @param q Logical. (starting with GDAL 1.7.0) To disable the progress bar on the console.
#' @param addalpha Logical. (starting with GDAL 1.7.0) Adds an alpha mask band to the VRT when the source raster have none. Mainly useful for RGB sources (or grey-level sources). The alpha band is filled on-the-fly with the value 0 in areas without any source raster, and with value 255 in areas with source raster. The effect is that a RGBA viewer will render the areas without source rasters as transparent and areas with source rasters as opaque. This option is not compatible with -separate.
#' @param hidenodata Logical. (starting with GDAL 1.7.0) Even if any band contains nodata value, giving this option makes the VRT band not report the NoData. Useful when you want to control the background color of the dataset. By using along with the -addalpha option, you can prepare a dataset which doesn't report nodata value but is transparent in areas with no data.
#' @param srcnodata Character. (starting with GDAL 1.7.0) Set nodata values for input bands (different values can be supplied for each band). If more than one value is supplied all values should be quoted to keep them together as a single operating system argument. If the option is not specified, the intrinsic nodata settings on the source datasets will be used (if they exist). The value set by this option is written in the NODATA element of each ComplexSource element. Use a value of None to ignore intrinsic nodata settings on the source datasets.
#' @param vrtnodata Character. (starting with GDAL 1.7.0) Set nodata values at the VRT band level (different values can be supplied for each band). If more than one value is supplied all values should be quoted to keep them together as a single operating system argument. If the option is not specified, intrinsic nodata settings on the first dataset will be used (if they exist). The value set by this option is written in the NoDataValue element of each VRTRasterBand element. Use a value of None to ignore intrinsic nodata settings on the source datasets.
#' @param a_srs Character. (starting with GDAL 1.10) Override the projection for the output file. The srs_def may be any of the usual GDAL/OGR forms, complete WKT, PROJ.4, EPSG:n or a file containing the WKT.
#' @param r Character. ("nearest" (default) | "bilinear" | "cubic" | "cubicspline" | "lanczos" | "average" | "mode"). (GDAL >= 2.0) Select a resampling algorithm.
#' @param input_file_list Character. To specify a text file with an input filename on each line.
#' @param overwrite Logical. Overwrite the VRT if it already exists.
## @param additional_commands Character. Additional commands to pass directly to ogrinfo.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @param ... Additional arguments.
#' 
#' @return NULL
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdalbuildvrt' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/gdalbuildvrt.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/gdalbuildvrt.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' layer1 <- system.file("external/tahoe_lidar_bareearth.tif", package="gdalUtils")
#' layer2 <- system.file("external/tahoe_lidar_highesthit.tif", package="gdalUtils")
#' output.vrt <- paste(tempfile(),".vrt",sep="")
#' gdalbuildvrt(gdalfile=c(layer1,layer2),output.vrt=output.vrt,separate=TRUE,verbose=TRUE)
#' gdalinfo(output.vrt)
#' }
#' @export

gdalbuildvrt <- function(gdalfile,output.vrt,
		tileindex,resolution,te,tr,tap,
		separate,b,sd,allow_projection_difference,q,
		addalpha,hidenodata,srcnodata,vrtnodata,
		a_srs,r,input_file_list,overwrite,
#		additional_commands,
		ignore.full_scan=TRUE,
		verbose=FALSE,
		...
)
{
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("tileindex","tap","addalpha",
							"hidenodata","separate","allow_projection_difference",
							"q","overwrite")),
			vector = list(
					varnames <- c("tr","te")),
			scalar = list(
					varnames <- c("b","sd")),
			character = list(
					varnames <- c("resolution","srcnodata","vrtnodata",
							"a_srs","r","input_file_list",
							"output.vrt","gdalfile")),
			repeatable = list(
					varnames <- NULL)
	)
	
	parameter_order <- c(
			"tileindex","resolution","te","tr","tap",
			"separate","b","sd","allow_projection_difference",
			"q","addalpha","hidenodata","srcnodata","vrtnodata",
			"a_srs","r","input_file_list","overwrite","output.vrt",
			"gdalfile")
	
	parameter_noflags <- c("output.vrt","gdalfile")
	
	parameter_doubledash <- NULL
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdalbuildvrt"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_noquotes=parameter_noquotes,
			parameter_doubledash=parameter_doubledash)
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	return(NULL)
}
