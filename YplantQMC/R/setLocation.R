#'Generate a location object
#'
#'@description Construct an object that contains a geographical location given by the
#'latitude, longitude, and (optionally), the longitude of the nearest timezone
#'border. This object is required in \code{\link{setMet}}.
#'
#'If the longitude of the nearest time zone border (\code{tzlong}) is not set,
#'it is assumed that all times are in 'local apparent time' (LAT), where
#'maximum solar altitude is reached at noon. If it is set, it is used to
#'calculate the time offset between clock time and solar time.
#'
#'If \code{lat} and \code{long} are not given, a location can be selected from
#'a simple map of the world. The \code{maps} package is needed for this
#'utility.
#'
#'You can also plot the location on a map of the world (see Examples).
#'
#'@param lat Latitude (degrees, southern hemisphere negative)
#'@param long Longitue (degrees)
#'@param tzlong Optional. Longitude of the nearest timezone border (degrees).
#'@return An object of class \code{yplocation}.
#'@author Remko Duursma
#'@seealso \code{\link{setMet}}
#'@keywords misc
#'@examples
#'
#'
#'# Set a location:
#'sydney <- setLocation(lat=-33.5, long=152, tzlong=150)
#'
#'\dontrun{
#'# Set a location from a map:
#'somewhere <- setLocation()
#'
#'# Plot locations ('maps' package required)
#'plot(sydney)
#'}
#'
#'@export
setLocation <- function(lat=NA, long=0, tzlong=NA){

	
	if(is.na(tzlong)){
		tzlong <- long
		tzlongset <- FALSE
	} else {
		tzlongset <- TRUE
	}
	
	if(is.na(lat)){
	  r <- require(maps)
	  if(!r) stop("Need to install package 'maps' to select location from a map.")
	  maps::map('world')
	  latlong <- locator(1)
	  lat <- latlong$y
	  long <- latlong$x
	  tzlong <- long
	  tzlongset <- FALSE
	  message("Time of sunrise will be local apparent time.")
	  message("Location successfully selected.")
	  flush.console()
	}
  
	l <- list()
	l$lat <- lat
	l$long <- long
	l$tzlong <- tzlong
	l$tzlongset <- tzlongset
	
	class(l) <- "yplocation"
	
return(l)	
}



