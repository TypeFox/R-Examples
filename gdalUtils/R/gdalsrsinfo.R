#' gdalsrsinfo
#' 
#' R wrapper for gdalsrsinfo: lists info about a given SRS in number of formats (WKT, PROJ.4, etc.)
#' 
#' @param srs_def Character. A raster dataset name. It can be either file name.
#' @param p Logical. Pretty-print where applicable (e.g. WKT).
#' @param V Logical. Validate SRS.
#' @param o Character. Output type ("default"|"all"|"wkt_all"|"proj4"|"wkt"|"wkt_simple"|"wkt_noct"|"wkt_esri"|"mapinfo"|"xml")
#' @param as.CRS Logical. Return a CRS object?  Default=FALSE.
## @param additional_commands Character. Additional commands to pass directly to gdalsrsinfo.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' 
#' @return character
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) and Matteo Mattiuzzi (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdalsrsinfo' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://www.gdal.org/gdalsrsinfo.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#' 
#' If as.CRS is set to TRUE, 'o' will automatically be set to "proj4" and the output
#' will be coerced to a CRS object for use with sp.
#'
#' @references \url{http://www.gdal.org/gdalsrsinfo.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' src_dataset <- system.file("external/tahoe_highrez.tif", package="gdalUtils")
#' # Command-line gdalsrsinfo call:
#' # gdalsrsinfo -o proj4 tahoe_highrez.tif
#' gdalsrsinfo(src_dataset,o="proj4",verbose=TRUE)
#' # Export as CRS:
#' gdalsrsinfo(src_dataset,as.CRS=TRUE,verbose=TRUE)
#' }
#' @import sp
#' @export

gdalsrsinfo <- function(srs_def,p,V,o,
#	additional_commands,
	as.CRS=FALSE,
	ignore.full_scan=TRUE,
	verbose=FALSE)
{
	if(as.CRS) o <- "proj4"

	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("p","V")),
			vector = list(
					varnames <- NULL),
			scalar = list(
					varnames <- c("sd")),
			character = list(
					varnames <- c("o","srs_def")),
			repeatable = list(
					varnames <- NULL)
	)
	
	parameter_order <- c(
			"p","V","o","srs_def")
	
	parameter_noflags <- c("srs_def")
	
	executable <- "gdalsrsinfo"
	# End gdalinfo setup
	
	cmd <- gdal_cmd_builder(
			executable=executable,
			parameter_variables=parameter_variables,
			parameter_values=parameter_values,
			parameter_order=parameter_order,
			parameter_noflags=parameter_noflags)
	
	if(verbose) message(paste("GDAL command being used:",cmd))
	
	cmd_output <- system(cmd,intern=TRUE) 
	
	if(verbose) { message(cmd_output) }
	if(as.CRS)
	{
		return(CRS(gsub("'","",cmd_output)))
	} else
	{
		return(cmd_output)
	}
}
