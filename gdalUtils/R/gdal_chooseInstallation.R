#' gdal_chooseInstallation
#' 
#' Choose a GDAL installation based on certain requirements.
#' 
#' @param hasDrivers Character. Which drivers must be available?
#' 
#' @return Numeric id of the most recent installation that matches the requirements. 
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net})
#' 
#' @details By default, the GDAL commands will use the installation found
#' at getOption("gdalUtils_gdalPath")[[1]], which is the most recent version
#' found on the system.  If the user has more than one GDAL installed (more
#' common on Windows and Mac systems than *nix systems), gdal_chooseInstallation
#' can be used to choose an installation (perhaps not the most recent one) that
#' has certain functionality, e.g. supports HDF4 formatted files.
#'
#' @references \url{http://www.gdal.org/gdal_translate.html}
#' @examples \dontrun{ 
#' # Choose the best installation that has both HDF4 and HDF5 drivers:
#' gdal_chooseInstallation(hasDrivers=c("HDF4","HDF5"))
#' # Get the version of this installation:
#' getOption("gdalUtils_gdalPath")[[
#' 	gdal_chooseInstallation(hasDrivers=c("HDF4","HDF5"))]]$version
#' }
#' @export

# TODO: hasWriteDrivers
# TODO: allow a user to force an ID
# TODO: minimum version requirement

gdal_chooseInstallation <- function(hasDrivers)
{	
	gdal_setInstallation()
	gdal_installations <- getOption("gdalUtils_gdalPath")
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	gdal_installation_ids <- seq(length(gdal_installations))
	current_match <- rep(x=TRUE,times=length(gdal_installations))
	
	if(!missing(hasDrivers))
	{
		driver_match <- sapply(gdal_installations,
				function(x,hasDrivers)
				{
					installation_drivers <- x$drivers
					if(!grepl("^2",x$version))
					{
						return(length(intersect(hasDrivers,installation_drivers$format_code)) 
										== length(hasDrivers))
					} else
					{
						# FIX FOR VERSION 2.0						
						matching_drivers <- grep(paste(hasDrivers,"-",sep=""),installation_drivers$format_code)
						return(length(matching_drivers)==length(hasDrivers))
					}				
				},
				hasDrivers=hasDrivers)	
		current_match <- current_match&driver_match
	}
	
	if(!any(current_match)) stop("No installations match.")
	
	installation_choice <- min(gdal_installation_ids[current_match])
	return(installation_choice)
}