#' gdaladdo
#' 
#' R wrapper for gdaladdo: builds or rebuilds overview images
#' 
#' @param filename Character. The file to build overviews for (or whose overviews must be removed).
#' @param levels Numeric. A list of integral overview levels to build. Ignored with clean=TRUE option.
#' @param r Character. ("nearest"|"average"|"gauss"|"cubic"|"average_mp"|"average_magphase"|"mode") Select a resampling algorithm.  Default is "nearest".
#' @param b Numeric. (available from GDAL 1.10) Select an input band band for overview generation. Band numbering starts from 1. Multiple -b switches may be used to select a set of input bands to generate overviews.
#' @param ro Logical. (available from GDAL 1.6.0) open the dataset in read-only mode, in order to generate external overview (for GeoTIFF especially).
#' @param clean Logical. (available from GDAL 1.7.0) remove all overviews.
#' @param oo Character. NAME=VALUE. (starting with GDAL 2.0) Dataset open option (format specific)
## @param additional_commands Character. Additional commands to pass directly to gdaladdo.
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  
#' 
#' @return NULL
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) (wrapper) and Frank Warmerdam (GDAL lead developer).
#' @details This is an R wrapper for the 'gdaladdo' function that is part of the 
#' Geospatial Data Abstraction Library (GDAL).  It follows the parameter naming
#' conventions of the original function, with some modifications to allow for more R-like
#' parameters.  For all parameters, the user can use a single character string following,
#' precisely, the gdalinfo format (\url{http://gdal.org/gdaladdo.html}), or,
#' in some cases, can use R vectors to achieve the same end.  
#' 
#' This function assumes the user has a working GDAL on their system.  If the 
#' "gdalUtils_gdalPath" option has been set (usually by gdal_setInstallation),
#' the GDAL found in that path will be used.  If nothing is found, gdal_setInstallation
#' will be executed to attempt to find a working GDAL.
#'
#' @references \url{http://www.gdal.org/gdaladdo.html}
#' 
#' @examples 
#' # We'll pre-check to make sure there is a valid GDAL install.
#' # Note this isn't strictly neccessary, as executing the function will
#' # force a search for a valid GDAL install.
#' gdal_setInstallation()
#' valid_install <- !is.null(getOption("gdalUtils_gdalPath"))
#' if(valid_install)
#' {
#' filename  <- system.file("external/tahoe_highrez.tif", package="gdalUtils")
#' temp_filename <- paste(tempfile(),".tif",sep="")
#' file.copy(from=filename,to=temp_filename,overwrite=TRUE)
#' gdalinfo(filename)
#' gdaladdo(r="average",temp_filename,levels=c(2,4,8,16),verbose=TRUE)
#' gdalinfo(temp_filename)
#' }
#' @export

gdaladdo <- function(filename,levels,
		r,b,ro,clean,oo,
#		additional_commands,
		ignore.full_scan=TRUE,
		verbose=FALSE)
{
	
	parameter_values <- as.list(environment())
	
	if(verbose) message("Checking gdal_installation...")
	gdal_setInstallation(ignore.full_scan=ignore.full_scan,verbose=verbose)
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	# Start gdalinfo setup
	parameter_variables <- list(
			logical = list(
					varnames <- c("ro","clean")),
			vector = list(
					varnames <- c("levels")),
			scalar = list(
					varnames <- c("")),
			character = list(
					varnames <- c("r","oo","filename")),
			repeatable = list(
					varnames <- c("b"))
	)
	
	parameter_order <- c(
			"r","b","ro","clean","oo","filename","levels")
	
	parameter_noflags <- c("filename","levels")
	
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
	
	return(NULL)
}
