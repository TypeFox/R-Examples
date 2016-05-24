#' gdal_rasterize
#' 
#' R wrapper for gdal_rasterize: burns vector geometries into a raster
#' 
#' @param src_datasource Character. Any OGR supported readable datasource.
#' @param dst_filename Character. The GDAL supported output file. Must support update mode access. Before GDAL 1.8.0, gdal_rasterize could not create new output files.
#' @param b Numeric. The band(s) to burn values into. Multiple -b arguments may be used to burn into a list of bands. The default is to burn into band 1.
#' @param i Logical. Invert rasterization. Burn the fixed burn value, or the burn value associated with the first feature into all parts of the image not inside the provided a polygon.
#' @param at Logical. Enables the ALL_TOUCHED rasterization option so that all pixels touched by lines or polygons will be updated not just those one the line render path, or whose center point is within the polygon. Defaults to disabled for normal rendering rules.
#' @param burn Numeric. A fixed value to burn into a band for all objects. A vector of burn options can be supplied, one per band being written to.
#' @param a Character. Identifies an attribute field on the features to be used for a burn in value. The value will be burned into all output bands.
#' @param threeD Logical. (GDAL parameter '3d') Indicates that a burn value should be extracted from the "Z" values of the feature. These values are adjusted by the burn value given by "-burn value" or "-a attribute_name" if provided. As of now, only points and lines are drawn in 3D. 
#' @param l Character. Indicates the layer(s) from the datasource that will be used for input features. May be specified multiple times, but at least one layer name or a -sql option must be specified.
#' @param where Character. An optional SQL WHERE style query expression to be applied to select features to burn in from the input layer(s).
#' @param sql Character. An SQL statement to be evaluated against the datasource to produce a virtual layer of features to be burned in.
#' @param dialect Character. (starting with GDAL 2.1) The SQL dialect. In some cases can be used to use (unoptimized) OGR SQL instead of the native SQL of an RDBMS by passing OGRSQL. Starting with GDAL 1.10, the "SQLITE" dialect can also be used with any datasource.
#' @param of Character. (GDAL >= 1.8.0) Select the output format. The default is GeoTIFF (GTiff). Use the short format name.
#' @param a_nodata Numeric. (GDAL >= 1.8.0) Assign a specified nodata value to output bands.
#' @param init Numeric. (GDAL >= 1.8.0) Pre-initialize the output image bands with these values. However, it is not marked as the nodata value in the output file. If only one value is given, the same value is used in all the bands.
#' @param a_srs Character. (GDAL >= 1.8.0) Override the projection for the output file. If not specified, the projection of the input vector file will be used if available. If incompatible projections between input and output files, no attempt will be made to reproject features. The srs_def may be any of the usual GDAL/OGR forms, complete WKT, PROJ.4, EPSG:n or a file containing the WKT.
#' @param co Character. (GDAL >= 1.8.0) Passes a creation option ("NAME=VALUE") to the output format driver. Multiple -co options may be listed. See format specific documentation for legal creation options for each format.
#' @param te Numeric. c(xmin,ymin,xmax,ymax) (GDAL >= 1.8.0) set georeferenced extents. The values must be expressed in georeferenced units. If not specified, the extent of the output file will be the extent of the vector layers.
#' @param tr Numeric. c(xres,yres) (GDAL >= 1.8.0) set target resolution. The values must be expressed in georeferenced units. Both must be positive values.
#' @param tap Logical. (GDAL >= 1.8.0) (target aligned pixels) align the coordinates of the extent of the output file to the values of the -tr, such that the aligned extent includes the minimum extent.
#' @param ts Numeric. c(width,height) (GDAL >= 1.8.0) set output file size in pixels and lines. Note that -ts cannot be used with -tr
#' @param ot Character. (GDAL >= 1.8.0) For the output bands to be of the indicated data type. Defaults to Float64
#' @param q Logical. (GDAL >= 1.8.0) Suppress progress monitor and other non-error output.
#' @param output_Raster Logical. Return output dst_filename as a RasterBrick?
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' @return NULL or if(output_Raster), a RasterBrick.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdal_rasterize' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalwarp format (\url{http://www.gdal.org/gdal_rasterize.html}), or,
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
#' @references \url{http://www.gdal.org/gdal_rasterize.html}
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(raster) && require(rgdal) && valid_install)
#' {
#' # Example from the original gdal_rasterize documentation:
#' # gdal_rasterize -b 1 -b 2 -b 3 -burn 255 -burn 0 
#' # 	-burn 0 -l tahoe_highrez_training tahoe_highrez_training.shp tempfile.tif
#' dst_filename_original  <- system.file("external/tahoe_highrez.tif", package="gdalUtils")
#' # Back up the file, since we are going to burn stuff into it.
#' dst_filename <- paste(tempfile(),".tif",sep="")
#' file.copy(dst_filename_original,dst_filename,overwrite=TRUE)
#' #Before plot:
#' plotRGB(brick(dst_filename))
#' src_dataset <- system.file("external/tahoe_highrez_training.shp", package="gdalUtils")
#' tahoe_burned <- gdal_rasterize(src_dataset,dst_filename,
#' 	b=c(1,2,3),burn=c(0,255,0),l="tahoe_highrez_training",verbose=TRUE,output_Raster=TRUE)
#' #After plot:
#' plotRGB(brick(dst_filename))
#' }
#' @import rgdal
#' @export

gdal_rasterize <- function(
		src_datasource,dst_filename,
		b,i,at,burn,a,threeD,l,where,sql,dialect,
		of,a_srs,co,a_nodata,init,
		te,tr,tap,ts,ot,q,
#		additional_commands,
		output_Raster=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE)
{
	if(output_Raster && (!requireNamespace("raster") || !requireNamespace("rgdal")))
	{
		warning("rgdal and/or raster not installed. Please install.packages(c('rgdal','raster')) or set output_Raster=FALSE")
		return(NULL)
	}
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Place all gdal function variables into these groupings:
	parameter_variables <- list(
			logical = list(
					varnames <- c(
							"i","at","threeD","tap","q"
					)),
			vector = list(
					varnames <- c(
						"init","te","tr","ts"	
					)),
			scalar = list(
					varnames <- c(
						"a_nodata"
					)),
			character = list(
					varnames <- c(
					"a","where","sql","dialect","of","a_srs","ot",
					"src_datasource","dst_filename"
					)),
			repeatable = list(
					varnames <- c(
						"b","burn","l","co"
					))
	)
	
	parameter_order <- c(
			"b","i","at","burn","a","threeD","l",
			"where","sql","dialect",
			"of","a_srs","co","a_nodata","init",
			"te","tr","tap","ts","ot","q",
			"src_datasource","dst_filename"
	)
	
	parameter_noflags <- c("src_datasource","dst_filename")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdal_rasterize"
	
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
		return(brick(dst_filename))	
	} else
	{
		return(NULL)
	}		
}