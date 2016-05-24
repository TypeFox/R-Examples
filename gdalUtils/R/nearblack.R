#' nearblack
#' 
#' R wrapper for nearblack: convert nearly black/white borders to black
#' 
#' @param infile Character. The input file. Any GDAL supported format, any number of bands, normally 8bit Byte bands.
#' @param o Character. outfile. The name of the output file to be created. Newly created files are created with the HFA driver by default (Erdas Imagine - .img)
#' @param of Character. format. (GDAL 1.8.0 or later) Select the output format. Use the short format name (GTiff for GeoTIFF for example).
#' @param co Character. "NAME=VALUE". (GDAL 1.8.0 or later) Passes a creation option to the output format driver. Multiple -co options may be listed. See format specific documentation for legal creation options for each format. Only valid when creating a new file.
#' @param white Logical. Search for nearly white (255) pixels instead of nearly black pixels.
#' @param color Numeric. c1,c2,c3...cn. (GDAL >= 1.9.0) Search for pixels near the specified color. May be specified multiple times. When -color is specified, the pixels that are considered as the collar are set to 0.
#' @param near Numeric. dist. Select how far from black, white or custom colors the pixel values can be and still considered near black, white or custom color. Defaults to 15.
#' @param nb Numeric. non_black_pixels. number of non-black pixels that can be encountered before the giving up search inwards. Defaults to 2.
#' @param setalpha Logical. (GDAL 1.8.0 or later) Adds an alpha band if the output file is specified and the input file has 3 bands, or sets the alpha band of the output file if it is specified and the input file has 4 bands, or sets the alpha band of the input file if it has 4 bands and no output file is specified. The alpha band is set to 0 in the image collar and to 255 elsewhere.
#' @param setmask Logical. (GDAL 1.8.0 or later) Adds a mask band to the output file, or adds a mask band to the input file if it does not already have one and no output file is specified. The mask band is set to 0 in the image collar and to 255 elsewhere.
#' @param q Logical. (GDAL 1.8.0 or later) Suppress progress monitor and other non-error output. 
## @param additional_commands Character. Additional commands to pass directly to gdaladdo.
#' @param output_Raster Logical. Return outfile as a RasterBrick?
#' @param overwrite Logical. If output file exists, OR if output file is not set (which would defualt to overwriting the input file), allow overwriting? 
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
## @param ... Other parameters to pass to nearblack.
#' 
#' @return NULL
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details #' This is an R wrapper for the 'nearblack' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/nearblack.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This utility will scan an image and try to set all pixels that are nearly or 
#' exactly black, white or one or more custom colors around the collar to black 
#' or white. This is often used to "fix up" lossy compressed airphotos so that 
#' color pixels can be treated as transparent when mosaicking.
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/nearblack.html}
#' 
#' @examples 
#'  # None available at present.
#' @export

nearblack <- function(
		infile,o,
		of,co,white,color,near,nb,setalpha,setmask,q,
#		additional_commands,
		output_Raster=FALSE,
		overwrite=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE#,
#		...
)
{
	
	if(output_Raster && (!requireNamespace("raster") || !requireNamespace("rgdal")))
	{
		warning("rgdal and/or raster not installed. Please install.packages(c('rgdal','raster')) or set output_Raster=FALSE")
		return(NULL)
	}
	
	if(missing(o) & !overwrite) stop("Warning: You are attempting to overwrite your input file. Set an output file 'o' or overwrite=T (to overwrite the input file) to proceed.")
	if(file.exists(o) & !overwrite) stop("Output file exists.  Set overwrite=T or pick another output name.")
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("white","setalpha","setmask","q")),
			vector = list(
					varnames <- c()),
			scalar = list(
					varnames <- c("near","nb")),
			character = list(
					varnames <- c("o","of","co","infile")),
			repeatable = list(
					varnames <- c("color"))
	)
	
	parameter_order <- c("white","setalpha","setmask","q",
			"near","nb",
			"o","of","co",
			"color",
			"infile"
	)
	
	parameter_noflags <- c("infile")
	
	parameter_doubledash <- NULL
	
	parameter_noquotes <- unlist(parameter_variables$vector)
	
	executable <- "gdaladdo"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags,
			parameter_doubledash=parameter_doubledash,
			parameter_noquotes=parameter_noquotes)
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(output_Raster)
	{
		if(missing(o)) 
		{ 
			return(brick(infile))
		} else
		{
			return(brick(o))
		}
	} else
	{
		return(NULL)
	}		
}
