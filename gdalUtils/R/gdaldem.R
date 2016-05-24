#' gdaldem
#' 
#' R wrapper for gdaldem: Tools to analyze and visualize DEMs. (since GDAL 1.7.0)
#' 
#' @param mode Character. ("hillshade"|"slope"|"aspect"|"color-relief"|"TRI"|"TPI"|"roughness")
#' @param input_dem Character. The input DEM raster to be processed.
#' @param output Character. The output raster produced.
#' @param of Character. Select the output format. The default is GeoTIFF (GTiff). Use the short format name.
#' @param compute_edges Logical. (GDAL >= 1.8.0) Do the computation at raster edges and near nodata values.
#' @param alg Character. "ZevenbergenThorne" (GDAL >= 1.8.0) Use Zevenbergen & Thorne formula, instead of Horn's formula, to compute slope & aspect. The litterature suggests Zevenbergen & Thorne to be more suited to smooth landscapes, whereas Horn's formula to perform better on rougher terrain.
#' @param b Numeric. Select an input band to be processed. Bands are numbered from 1.
#' @param co Character. (GDAL >= 1.8.0) Passes a creation option ("NAME=VALUE") to the output format driver. Multiple -co options may be listed. See format specific documentation for legal creation options for each format.
#' @param q Logical. Suppress progress monitor and other non-error output.
#' @param z Numeric. (mode=="hillshade") vertical exaggeration used to pre-multiply the elevations.
#' @param s Numeric. (mode=="hillshade" | mode=="slope) ratio of vertical units to horizontal. If the horizontal unit of the source DEM is degrees (e.g Lat/Long WGS84 projection), you can use scale=111120 if the vertical units are meters (or scale=370400 if they are in feet).
#' @param az Numeric. (mode=="hillshade") azimuth of the light, in degrees. 0 if it comes from the top of the raster, 90 from the east, ... The default value, 315, should rarely be changed as it is the value generally used to generate shaded maps.
#' @param alt Numeric. (mode=="hillshade") altitude of the light, in degrees. 90 if the light comes from above the DEM, 0 if it is raking light.
#' @param combined Character. (mode=="hillshade") "combined shading" (starting with GDAL 1.10) a combination of slope and oblique shading.
#' @param p Logical. (mode=="slope") if specified, the slope will be expressed as percent slope. Otherwise, it is expressed as degrees.
#' @param trigonometric Logical. (mode=="aspect") return trigonometric angle instead of azimuth. Thus 0deg means East, 90deg North, 180deg West, 270deg South.
#' @param zero_for_flat Logical. (mode=="aspect") By using those 2 options, the aspect returned by gdaldem aspect should be identical to the one of GRASS r.slope.aspect. Otherwise, it's identical to the one of Matthew Perry's aspect.cpp utility.
#' @param color_text_file Character. (mode=="color-relief") text-based color configuration file (see Description).
#' @param alpha Logical. (mode=="color-relief") add an alpha channel to the output raster.
#' @param exact_color_entry Logical. (mode=="color-relief") use strict matching when searching in the color configuration file. If none matching color entry is found, the "0,0,0,0" RGBA quadruplet will be used.
#' @param nearest_color_entry Logical. (mode=="color-relief") use the RGBA quadruplet corresponding to the closest entry in the color configuration file.
## @param additional_commands Character. Additional commands to pass directly to gdalsrsinfo.
#' @param output_Raster Logical. Return output dst_dataset as a RasterBrick?
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' 
#' @return NULL or if(output_Raster), a RasterBrick.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Matthew Perry, Even Rouault, Howard Butler, and Chris Yesson  (GDAL developers).
#' @details This is an R wrapper for the 'gdaldem' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://www.gdal.org/gdaldem.html}), or,
#' in some cases, use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#' 
#' The user can choose to (optionally) return a RasterBrick of the output file (assuming
#' raster/rgdal supports the particular output format).
#'
#' @references \url{http://www.gdal.org/gdaldem.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install
#' # and that raster and rgdal are also installed.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(require(raster) && require(rgdal) && valid_install)
#' {
#' # We'll pre-check for a proper GDAL installation before running these examples:
#' gdal_setInstallation()
#' if(!is.null(getOption("gdalUtils_gdalPath")))
#' {
#' input_dem  <- system.file("external/tahoe_lidar_highesthit.tif", package="gdalUtils")
#' plot(raster(input_dem),col=gray.colors(256))
#' 
#' # Hillshading:
#' # Command-line gdaldem call:
#' # gdaldem hillshade tahoe_lidar_highesthit.tif output_hillshade.tif
#' output_hillshade <- gdaldem(mode="hillshade",input_dem=input_dem,
#'	output="output_hillshade.tif",output_Raster=TRUE,verbose=TRUE)
#' plot(output_hillshade,col=gray.colors(256))
#' 
#' # Slope:
#' # Command-line gdaldem call:
#' # gdaldem slope tahoe_lidar_highesthit.tif output_slope.tif -p
#' output_slope <- gdaldem(mode="slope",input_dem=input_dem,
#'	output="output_slope.tif",p=TRUE,output_Raster=TRUE,verbose=TRUE)
#' plot(output_slope,col=gray.colors(256))
#' 
#' # Aspect:
#' # Command-line gdaldem call:
#' # gdaldem aspect tahoe_lidar_highesthit.tif output_aspect.tif
#' output_aspect <- gdaldem(mode="aspect",input_dem=input_dem,
#'	output="output_aspect.tif",output_Raster=TRUE,verbose=TRUE)
#' plot(output_aspect,col=gray.colors(256))
#' }
#' }
#' @import rgdal
#' @export

# TODO: Fully document this.

gdaldem <- function(
		mode,input_dem,output,
		of,compute_edges,alg,b,co,q,
		z,s,az,alt,combined,
		p,
		trigonometric,zero_for_flat,
		color_text_file,alpha,exact_color_entry,nearest_color_entry,
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
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("compute_edges","q",
							"p",
							"trigonometric","zero_for_flat",
							"alpha","exact_color_entry","nearest_color_entry")),
			vector = list(
					varnames <- c("b")),
			scalar = list(
					varnames <- c("z","s","az","alt")),
			character = list(
					varnames <- c("mode","input_dem","output","of","alg",
							"combined","color_text_file")),
			repeatable = list(
					varnames <- c("co"))
	)
	
	parameter_order <- c("mode",
			"input_dem","color_text_file","output",
			"of","compute_edges","alg","b","co","q",
			"z","s","az","alt","combined",
			"p",
			"trigonometric","zero_for_flat",
			"alpha","exact_color_entry","nearest_color_entry"
			)
	
	parameter_noflags <- c("mode","input_dem","color_text_file","output")
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdaldem"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_noquotes=parameter_noquotes)
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(verbose) { message(cmd_output) }
	if(output_Raster)
	{
		return(brick(output))	
	} else
	{
		return(NULL)
	}		
}
