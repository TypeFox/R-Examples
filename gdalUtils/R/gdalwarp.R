#' gdalwarp
#' 
#' R wrapper for gdalwarp: image reprojection and warping utility
#' 
#' @param srcfile Character. The source file name(s).
#' @param dstfile Character. The destination file name.
#' @param s_srs Character. source spatial reference set. The coordinate systems that can be passed are anything supported by the OGRSpatialReference.SetFromUserInput() call, which includes EPSG PCS and GCSes (ie. EPSG:4296), PROJ.4 declarations (as above), or the name of a .prf file containing well known text.
#' @param t_srs Character. target spatial reference set. The coordinate systems that can be passed are anything supported by the OGRSpatialReference.SetFromUserInput() call, which includes EPSG PCS and GCSes (ie. EPSG:4296), PROJ.4 declarations (as above), or the name of a .prf file containing well known text.
#' @param to Character. set a transformer option suitable to pass to GDALCreateGenImgProjTransformer2().
#' @param order Numeric. order of polynomial used for warping (1 to 3). The default is to select a polynomial order based on the number of GCPs.
#' @param tps Logical. Force use of thin plate spline transformer based on available GCPs.
#' @param rpc Logical. Force use of RPCs. 
#' @param geoloc Logical. Force use of Geolocation Arrays.
#' @param et Numeric. error threshold for transformation approximation (in pixel units - defaults to 0.125).
#' @param refine_gcps Numeric. (GDAL >= 1.9.0) refines the GCPs by automatically eliminating outliers. Outliers will be eliminated until minimum_gcps are left or when no outliers can be detected. The tolerance is passed to adjust when a GCP will be eliminated. Note that GCP refinement only works with polynomial interpolation. The tolerance is in pixel units if no projection is available, otherwise it is in SRS units. If minimum_gcps is not provided, the minimum GCPs according to the polynomial model is used.
#' @param te Numeric. (c(xmin,ymin,xmax,ymax)). set georeferenced extents of output file to be created (in target SRS).
#' @param te_srs Character. srs_def. (GDAL >= 2.0) Specifies the SRS in which to interpret the coordinates given with -te. The srs_def may be any of the usual GDAL/OGR forms, complete WKT, PROJ.4, EPSG:n or a file containing the WKT. This must not be confused with -t_srs which is the target SRS of the output dataset. -te_srs is a conveniency e.g. when knowing the output coordinates in a geodetic long/lat SRS, but still wanting a result in a projected coordinate system.
#' @param tr Numeric. (c(xres,yres)). set output file resolution (in target georeferenced units)
#' @param tap Logical. (GDAL >= 1.8.0) (target aligned pixels) align the coordinates of the extent of the output file to the values of the -tr, such that the aligned extent includes the minimum extent.
#' @param ts Numeric. (c(width,height)). set output file size in pixels and lines. If width or height is set to 0, the other dimension will be guessed from the computed resolution. Note that -ts cannot be used with -tr
#' @param ovr Character. (level | "AUTO" | "AUTO-n" | "NONE"). (GDAL >= 2.0) To specify which overview level of source files must be used. The default choice, AUTO, will select the overview level whose resolution is the closest to the target resolution. Specify an integer value (0-based, i.e. 0=1st overview level) to select a particular level. Specify AUTO-n where n is an integer greater or equal to 1, to select an overview level below the AUTO one. Or specify NONE to force the base resolution to be used.
#' @param wo Character. Set a warp options. The GDALWarpOptions::papszWarpOptions docs show all options. Multiple -wo options may be listed.
#' @param ot Character. For the output bands to be of the indicated data type.
#' @param wt Character. Working pixel data type. The data type of pixels in the source image and destination image buffers.
#' @param r Character. resampling_method. ("near"|"bilinear"|"cubic"|"cubicspline"|"lanczos"|"average"|"mode"|"max"|"min"|"med"|"q1"|"q3")  See Description.
#' @param srcnodata Character. Set nodata masking values for input bands (different values can be supplied for each band). If more than one value is supplied all values should be quoted to keep them together as a single operating system argument. Masked values will not be used in interpolation. Use a value of None to ignore intrinsic nodata settings on the source dataset.
#' @param dstnodata Character. Set nodata values for output bands (different values can be supplied for each band). If more than one value is supplied all values should be quoted to keep them together as a single operating system argument. New files will be initialized to this value and if possible the nodata value will be recorded in the output file. Use a value of None to ensure that nodata is not defined (GDAL>=2.0). If this argument is not used then nodata values will be copied from the source dataset (GDAL>=2.0).
#' @param dstalpha Logical. Create an output alpha band to identify nodata (unset/transparent) pixels.
#' @param wm Numeric. Set the amount of memory (in megabytes) that the warp API is allowed to use for caching. 
#' @param multi Logical. Use multithreaded warping implementation. Multiple threads will be used to process chunks of image and perform input/output operation simultaneously.
#' @param q Logical. Be quiet.
#' @param of Character. Select the output format. The default is GeoTIFF (GTiff). Use the short format name.
#' @param co Character. passes a creation option to the output format driver. Multiple -co options may be listed. See format specific documentation for legal creation options for each format.
#' @param cutline Character. Enable use of a blend cutline from the name OGR support datasource.
#' @param cl Character. Select the named layer from the cutline datasource.
#' @param cwhere Character. Restrict desired cutline features based on attribute query.
#' @param csql Character. Select cutline features using an SQL query instead of from a layer with -cl.
#' @param cblend Numeric. Set a blend distance to use to blend over cutlines (in pixels).
#' @param crop_to_cutline Logical. (GDAL >= 1.8.0) Crop the extent of the target dataset to the extent of the cutline.
#' @param overwrite Logical. (GDAL >= 1.8.0) Overwrite the target dataset if it already exists.
#' @param nomd Logical. (GDAL >= 1.10.0) Do not copy metadata. Without this option, dataset and band metadata (as well as some band information) will be copied from the first source dataset. Items that differ between source datasets will be set to * (see -cvmd option).
#' @param cvmd Character. (GDAL >= 1.10.0) Value to set metadata items that conflict between source datasets (default is "*"). Use "" to remove conflicting items.
#' @param setci Logical. (GDAL >= 1.10.0) Set the color interpretation of the bands of the target dataset from the source dataset.
#' @param oo Character. NAME=VALUE. (starting with GDAL 2.0) Dataset open option (format specific).
#' @param doo Character. NAME=VALUE. (starting with GDAL 2.1) Output dataset open option (format specific).
## @param additional_commands Character. Additional commands to pass directly to gdalwarp.
#' @param output_Raster Logical. Return output dst_dataset as a RasterBrick?
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @param ... Additional arguments.
#' 
#' @return NULL or if(output_Raster), a RasterBrick.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdalwarp' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalwarp format (\url{http://www.gdal.org/gdalwarp.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL that has the right drivers 
#' as specified with the "of" (output format) parameter.
#' 
#' The resampling_methods available are as follows:
#' \itemize{
#' \item{near: nearest neighbour resampling (default, fastest algorithm, worst interpolation quality).}
#' \item{bilinear: bilinear resampling.}
#' \item{cubic: cubic resampling.}
#' \item{cubicspline: cubic spline resampling.}
#' \item{lanczos: Lanczos windowed sinc resampling.}
#' \item{average: average resampling, computes the average of all non-NODATA contributing pixels. (GDAL >= 1.10.0)}
#' \item{mode: mode resampling, selects the value which appears most often of all the sampled points. (GDAL >= 1.10.0)}
#' \item{max: maximum resampling, selects the maximum value from all non-NODATA contributing pixels. (GDAL >= 2.0.0)}
#' \item{min: minimum resampling, selects the minimum value from all non-NODATA contributing pixels. (GDAL >= 2.0.0)}
#' \item{med: median resampling, selects the median value of all non-NODATA contributing pixels. (GDAL >= 2.0.0)}
#' \item{q1: first quartile resampling, selects the first quartile value of all non-NODATA contributing pixels. (GDAL >= 2.0.0)}
#' \item{q3: third quartile resampling, selects the third quartile value of all non-NODATA contributing pixels. (GDAL >= 2.0.0)}
#' }

#' The user can choose to (optionally) return a RasterBrick of the output file (assuming
#' raster/rgdal supports the particular output format).
#'
#' @references \url{http://www.gdal.org/gdalwarp.html}
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(raster) && require(rgdal) && valid_install)
#' {
#' # Example from the original gdal_translate documentation:
#' src_dataset <- system.file("external/tahoe_highrez.tif", package="gdalUtils")
#' # Command-line gdalwarp call:
#' # gdalwarp -t_srs '+proj=utm +zone=11 +datum=WGS84' raw_spot.tif utm11.tif
#' gdalwarp(src_dataset,dstfile="tahoe_highrez_utm11.tif",
#' 		t_srs='+proj=utm +zone=11 +datum=WGS84',output_Raster=TRUE,
#' 		overwrite=TRUE,verbose=TRUE)
#' }
#' @import rgdal
#' @export

gdalwarp <- function(
		#help_general,formats, # Need to fix these
		srcfile,dstfile,
		s_srs,t_srs,to,
		order,tps,rpc,geoloc,et,refine_gcps,te,te_srs,tr,tap,ts,ovr,wo,ot,wt,r,srcnodata,dstnodata,
		dstalpha,wm,multi,q,of="GTiff",co,cutline,cl,cwhere,csql,cblend,crop_to_cutline,
		overwrite,nomd,cvmd,setci,oo,doo,
#		additional_commands,
		output_Raster=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE,
		...)
{
	if(output_Raster && (!requireNamespace("raster") || !requireNamespace("rgdal")))
	{
		warning("rgdal and/or raster not installed. Please install.packages(c('rgdal','raster')) or set output_Raster=FALSE")
		return(NULL)
	}
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Place all gdal function variables into these groupings:
	parameter_variables <- list(
			logical = list(
					varnames <- c(
							"tps","rpc","geoloc","tap","dstalpha",
							"multi","q","crop_to_cutline","overwrite","nomd",
							"setci"
					)),
			vector = list(
					varnames <- c(
							"te","tr","ts"
					)),
			scalar = list(
					varnames <- c(
							"order","et","refine_gcps","wm",
							"cblend"
					)),
			character = list(
					varnames <- c(
							"s_srs","t_srs","to","te_srs","ovr","ot","wt","r",
							"srcnodata","dstnodata","of","cutline","cl",
							"cwhere","csql","cvmd","oo","doo","dstfile"
					)),
			repeatable = list(
					varnames <- c(
							"wo","co","srcfile"
					))
	)
	
	parameter_order <- c(
			"tps","rpc","geoloc","tap","dstalpha",
			"multi","q","crop_to_cutline","overwrite","nomd",
			"setci",
			"te","te_srs","tr","ts","ovr",
			"order","et","refine_gcps","wm",
			"cblend",
			"s_srs","t_srs","to","ot","wt","r",
			"srcnodata","dstnodata","of","cutline","cl",
			"cwhere","csql","cvmd",
			"wo","co","oo","doo",
			"srcfile","dstfile"
	)
	
	parameter_noflags <- c("srcfile","dstfile")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdalwarp"
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_noquotes=parameter_noquotes,
			gdal_installation_id=gdal_chooseInstallation(hasDrivers=of))
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(output_Raster)
	{
		return(brick(dstfile))	
	} else
	{
		return(NULL)
	}		
}