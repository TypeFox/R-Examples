## ###########################################################
## PURPOSE: Project coordinates from longitude/latitude to model units.
##
## INPUTS:
##   longitude: vector of longitudes for the points to be projected
##   latitude: vector of latitudes for the points to be projected
##   file: File name of Models3-formatted file giving the desired model
##     projection.
##   units: Units to be used for the projected coordinates; that is,
##     either "m" (meters) or "km" (kilometers). The option "deg" does
##     not make sense because we would not need to project coordinates
##     from a long/lat reference system to a Models3-formatted file
##     which was also gridded by long/lat.  If unspecified, the default
##     is "km".
##   ...: Other arguments to pass to get.proj.info.M3 function.
##     In this case, the only relevant argument would be the earth
##     radius to use when doing the projections.
##
## RETURNS: A list containing two elements "coords" and "units".  The
##   element "coords" contains a matrix with the projected coordinates.
##   The element "units" contains the units of the coordinates ("km" or "m").
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-05-19
## ###########################################################
project.lonlat.to.M3 <- function(longitude, latitude, file,
                                      units, ...){

  ## Form projection string describing the projection in the given
  ## Models3-formatted file.
  proj.string <- get.proj.info.M3(file, ...)


  ## Return an error if the Models3-formatted file given is already
  ## referenced by long/lat.  If it is, then we don't need to project
  ## the given coordinates, since they're already on the long/lat system.
  if ( substring(proj.string, first=7, last=13)=="longlat" )
    stop(paste("No need to project, since file ", file,
               " is gridded on long/lat system.", sep=""))


  ## Project locations from longitude and latitude onto CMAQ units (by
  ## default, this is done in meters).
  coords.proj <- project(cbind(longitude, latitude), proj=proj.string)
  colnames(coords.proj) <- c("x", "y")


  ## ##########################
  ## If the user does not specify desired units, then use "km".  If
  ## user specifies an option other than "km" or "m", give message and
  ## exit function.

  if (missing(units)){
    coords.proj <- coords.proj/1000
    units <- "km"
    }
  else if (units=="km"){
    coords.proj <- coords.proj/1000
    units <- "km"
  }
  else if (units=="m")
    units <- "m"
  else
    stop(paste(units, " is not a valid option for 'units'.", sep=""))
  ## ##########################


  ## Return a list, with the coords in the first position and the
  ## units of those coords in the second.
  x <- list(coords=coords.proj, units=units)
  return(x)
}
