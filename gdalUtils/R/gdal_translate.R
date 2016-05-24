#' gdal_translate
#' 
#' R wrapper for gdal_translate
#' 
#' @param src_dataset Character. The source dataset name. It can be either file name, URL of data source or subdataset name for multi-dataset files.
#' @param dst_dataset Character. The destination file name.
#' @param ot Character. ("Byte"/"Int16"/"UInt16"/"UInt32"/"Int32"/"Float32"/"Float64"/"CInt16"/"CInt32"/"CFloat32"/"CFloat64"). For the output bands to be of the indicated data type.
#' @param strict Logical. Don't be forgiving of mismatches and lost data when translating to the output format.
#' @param of Character. Select the output format. The default is GeoTIFF (GTiff). Use the short format name.
#' @param b Numeric or Character. Select an input band band for output. Bands are numbered from 1. Multiple bands may be used to select a set of input bands to write to the output file, or to reorder bands. Starting with GDAL 1.8.0, band can also be set to "mask,1" (or just "mask") to mean the mask band of the first band of the input dataset.
#' @param mask Numeric. (GDAL >= 1.8.0) Select an input band band to create output dataset mask band. Bands are numbered from 1. band can be set to "none" to avoid copying the global mask of the input dataset if it exists. Otherwise it is copied by default ("auto"), unless the mask is an alpha channel, or if it is explicitly used to be a regular band of the output dataset ("-b mask"). band can also be set to "mask,1" (or just "mask") to mean the mask band of the 1st band of the input dataset.
#' @param expand Character. ("gray"|"rgb"|"rgba").  (From GDAL 1.6.0) To expose a dataset with 1 band with a color table as a dataset with 3 (RGB) or 4 (RGBA) bands. Useful for output drivers such as JPEG, JPEG2000, MrSID, ECW that don't support color indexed datasets. The 'gray' value (from GDAL 1.7.0) enables to expand a dataset with a color table that only contains gray levels to a gray indexed dataset.
#' @param outsize Numeric. (c(xsize[percentage],ysize[percentage])). Set the size of the output file. Outsize is in pixels and lines unless '\%' is attached in which case it is as a fraction of the input image size.
#' @param tr Numeric. c(xres,yres). (starting with GDAL 2.0) set target resolution. The values must be expressed in georeferenced units. Both must be positive values. This is exclusive with -outsize and -a_ullr.
#' @param r Character. resampling_method. ("nearest"|"bilinear"|"cubic"|"cubicspline"|"lanczos"|"average"|"mode")  (GDAL >= 2.0) Select a resampling algorithm.
#' @param scale Numeric. (c(src_min,src_max,dst_min,dst_max)). Rescale the input pixels values from the range src_min to src_max to the range dst_min to dst_max. If omitted the output range is 0 to 255. If omitted the input range is automatically computed from the source data.
#' @param exponent Numeric. (From GDAL 1.11) To apply non-linear scaling with a power function. exp_val is the exponent of the power function (must be postive). This option must be used with the -scale option. If specified only once, -exponent applies to all bands of the output image. It can be repeated several times so as to specify per band parameters. It is also possible to use the "-exponent_bn" syntax where bn is a band number (e.g. "-exponent_2" for the 2nd band of the output dataset) to specify the parameters of one or several specific bands.
#' @param unscale Logical. Apply the scale/offset metadata for the bands to convert scaled values to unscaled values. It is also often necessary to reset the output datatype with the -ot switch.
#' @param srcwin Numeric. (c(xoff,yoff,xsize,ysize)).  Selects a subwindow from the source image for copying based on pixel/line location.
#' @param projwin Numeric. (c(ulx,uly,lrx,lry)).  Selects a subwindow from the source image for copying (like -srcwin) but with the corners given in georeferenced coordinates.
#' @param projwin_srs Character. srs_def. (GDAL >= 2.0) Specifies the SRS in which to interpret the coordinates given with -projwin. The srs_def may be any of the usual GDAL/OGR forms, complete WKT, PROJ.4, EPSG:n or a file containing the WKT. Note that this does not cause reprojection of the dataset to the specified SRS.
#' @param epo Logical. (Error when Partially Outside)  (GDAL >= 1.10) If this option is set, -srcwin or -projwin values that falls partially outside the source raster extent will be considered as an error. The default behaviour starting with GDAL 1.10 is to accept such requests, when they were considered as an error before.
#' @param eco Logical. (Error when Completely Outside) (GDAL >= 1.10) Same as -epo, except that the criterion for erroring out is when the request falls completely outside the source raster extent.
#' @param a_srs Character.  Override the projection for the output file. The srs_def may be any of the usual GDAL/OGR forms, complete WKT, PROJ.4, EPSG:n or a file containing the WKT.
#' @param a_ullr Numeric. (c(ulx,uly,lrx,lry)). Assign/override the georeferenced bounds of the output file. This assigns georeferenced bounds to the output file, ignoring what would have been derived from the source file.
#' @param a_nodata Numeric. Assign a specified nodata value to output bands. Starting with GDAL 1.8.0, can be set to none to avoid setting a nodata value to the output file if one exists for the source file
#' @param mo Character. ("META-TAG=VALUE").  Passes a metadata key and value to set on the output dataset if possible.
#' @param co Character. ("NAME=VALUE"). Passes a creation option to the output format driver. Multiple -co options may be listed. See format specific documentation for legal creation options for each format.
#' @param gcp Numeric. (c(pixel,line,easting,northing(,elevation))). Add the indicated ground control point to the output dataset. This option may be provided multiple times to provide a set of GCPs.
#' @param q Logical. Suppress progress monitor and other non-error output.
#' @param sds Logical. Copy all subdatasets of this file to individual output files. Use with formats like HDF or OGDI that have subdatasets.
#' @param stats Logical. (GDAL >= 1.8.0) Force (re)computation of statistics.
#' @param norat Logical. (GDAL >= 1.11) Do not copy source RAT into destination dataset.
#' @param oo Character. NAME=VALUE. (starting with GDAL 2.0) Dataset open option (format specific)
#' @param sd_index Numeric. If the file is an HDF4 or HDF5 file, which subdataset should be returned (1 to the number of subdatasets)?  If this flag is used, src_dataset should be the filename of the multipart file.  This parameter only works if the subdataset names follow the SUBDATASET_n_NAME convention.
#' @param output_Raster Logical. Return output dst_dataset as a RasterBrick?
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @param ... Additional arguments.
#' 
#' @return NULL or if(output_Raster), a RasterBrick.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdal_translate' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdal_translate format (\url{http://www.gdal.org/gdal_translate.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL that has the right drivers 
#' as specified with the "of" (output format) parameter.
#' 
#' The user can choose to (optionally) return a RasterBrick of the output file (assuming
#' raster/rgdal supports the particular output format).
#'
#' @references \url{http://www.gdal.org/gdal_translate.html}
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
#' # Original gdal_translate call:
#' # gdal_translate -of GTiff -co "TILED=YES" tahoe_highrez.tif tahoe_highrez_tiled.tif
#' gdal_translate(src_dataset,"tahoe_highrez_tiled.tif",of="GTiff",co="TILED=YES",verbose=TRUE)
#' # Pull out a chunk and return as a raster:
#' gdal_translate(src_dataset,"tahoe_highrez_tiled.tif",of="GTiff",
#' srcwin=c(1,1,100,100),output_Raster=TRUE,verbose=TRUE)
#' # Notice this is the equivalent, but follows gdal_translate's parameter format:
#' gdal_translate(src_dataset,"tahoe_highrez_tiled.tif",of="GTiff",
#' srcwin="1 1 100 100",output_Raster=TRUE,verbose=TRUE)
#' }
#' \dontrun{ 
#' # Extract the first subdataset from an HDF4 file:
#' hdf4_dataset <- system.file("external/test_modis.hdf", package="gdalUtils")
#' gdal_translate(hdf4_dataset,"test_modis_sd1.tif",sd_index=1)
#' }
#' @import rgdal
#' @export

# TODO: return all subdatasets if sds=TRUE as a list of bricks
# TODO: add in overwrite=TRUE/FALSE capabilties in this and other functions.
#		(Right now these are fairly unsafe, as they can easily overwrite and
#		existing file).
# TODO: return > 1 sd

gdal_translate <- function(src_dataset,dst_dataset,ot,strict,of="GTiff",
		b,mask,expand,outsize,tr,r,scale,exponent,unscale,srcwin,projwin,
		projwin_srs,epo,eco,
		a_srs,a_ullr,a_nodata,mo,co,gcp,q,sds,stats,norat,oo,
#		additional_commands,
		sd_index,
		output_Raster=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE,
		...
)
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
	
	if(!missing(sd_index))
	{
		parameter_values$src_dataset <- get_subdatasets(src_dataset,names_only=TRUE)[sd_index]
	}
	
	parameter_variables <- list(
			logical = list(
					varnames <- c("strict","unscale","epo","eco","q","sds","stats","norat")),
			vector = list(
					varnames <- c("outsize","tr","scale","exponent","srcwin","projwin","a_ullr","gcp")),
			scalar = list(
					varnames <- c("a_nodata")),
			character = list(
					varnames <- c("ot","of","mask","expand","r","projwin_srs","a_srs","oo","src_dataset","dst_dataset")),
			repeatable = list(
					varnames <- c("b","mo","co")))
	
	parameter_order <- c(
			"strict","exponent","unscale","epo","eco","q","sds","stats",
			"norat",
			"outsize","tr","scale","srcwin","projwin","a_ullr","gcp",
			"a_nodata",
			"ot","of","mask","expand","r","projwin_srs","a_srs",
			"b","mo","co","oo",
			"src_dataset","dst_dataset")
	
	parameter_noflags <- c("src_dataset","dst_dataset")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
		
	executable <- "gdal_translate"
	
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
	
	if(verbose) { message(cmd_output) } 
	
	# (Optional) return Raster
	if(output_Raster)
	{
		if(missing(sds)) 
		{
			return(brick(dst_dataset))
		}
		else
		{
			if(!sds) 
			{
				return(brick(dst_dataset))
			}
			else 
			{
				return(NULL)
			}
		}
	}
	else(return(NULL))
}
