#' gdal_contour
#' 
#' R wrapper for gdal_contour: builds vector contour lines from a raster elevation model
#' 
#' @param src_filename Character. Any OGR supported readable datasource.
#' @param dst_filename Character. The OGR supported output file. 
#' @param b Numeric. Picks a particular band to get the DEM from. Defaults to band 1.
#' @param a Character. Provides a name for the attribute in which to put the elevation. If not provided no elevation attribute is attached.
#' @param threeD Logical. (GDAL parameter '3d') Force production of 3D vectors instead of 2D. Includes elevation at every vertex.
#' @param inodata Logical. Ignore any nodata value implied in the dataset - treat all values as valid.
#' @param snodata Numeric. Input pixel value to treat as "nodata".
#' @param f Character. Create output in a particular format, default is "ESRI Shapefiles".
#' @param dsco Character. Dataset creation option (format specific).  Follows "NAME=VALUE" format.
#' @param lco Character. Layer creation option (format specific).  Follows "NAME=VALUE" format.
#' @param i Numeric. Elevation interval between contours.
#' @param off Numeric. Offset from zero relative to which to interpret intervals.
#' @param fl Character. Name one or more "fixed levels" to extract.
#' @param nln Character. Provide a name for the output vector layer. Defaults to "contour". 
#' @param output_Vector Logical. Return output dst_filename as a Spatial* object.  Currently only works with f="ESRI Shapefile".
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  

#' @return output vector filename or SpatialLinesDataFrame object.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdal_contour' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdal_contour format (\url{http://www.gdal.org/gdal_contour.html}), or,
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
#' @references \url{http://www.gdal.org/gdal_contour.html}
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(raster) && require(rgdal) && valid_install)
#' {
#' # Example from the original gdal_contour documentation:
#' # 	gdal_contour -a elev dem.tif contour.shp -i 10.0 
#' # Choose a DEM:
#' input_dem  <- system.file("external/tahoe_lidar_bareearth.tif", package="gdalUtils")
#' # Setup an output filename (shapefile):
#' output_shapefile <- paste(tempfile(),".shp",sep="")
#' contour_output <- gdal_contour(src_filename=input_dem,dst_filename=output_shapefile,
#' 		a="Elevation",i=5.,output_Vector=TRUE)
#' # Plot the contours using spplot:
#' spplot(contour_output["Elevation"],contour=TRUE)
#' }
#' @import rgdal
#' @export

gdal_contour <- function(
		src_filename,dst_filename,
		b,a,threeD,inodata,snodata,i,
		f="ESRI Shapefile",
		dsco,lco,off,fl,nln,
		output_Vector=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE)
{
	if(output_Vector && (!requireNamespace("rgdal")))
	{
		warning("rgdal not installed. Please install.packages('rgdal') or set output_Vector=FALSE")
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
							"threeD","inodata"
					)),
			vector = list(
					varnames <- c(
					
					)),
			scalar = list(
					varnames <- c(
							"b","snodata","i","off"
					)),
			character = list(
					varnames <- c(
							"a","f","fl","nln","src_filename","dst_filename"
					)),
			repeatable = list(
					varnames <- c(
							"dsco","lco"
					))
	)
	
	parameter_order <- c(
			"b","a","threeD","inodata","snodata",
			"f","dsco","lco","i","off","fl","nln",
			"src_filename","dst_filename"
	)
	
	parameter_noflags <- c("src_filename","dst_filename")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdal_contour"
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_noquotes=parameter_noquotes,
			#		gdal_installation_id=gdal_chooseInstallation(hasDrivers=of))
			gdal_installation_id=gdal_chooseInstallation())
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(output_Vector && f=="ESRI Shapefile")
	{
		return(
				readOGR(dsn=dirname(dst_filename),
						layer=remove_file_extension(basename(dst_filename)))
		)
		#	return(brick(dst_filename))	
	} else
	{
		return(dst_filename)
	}		
}