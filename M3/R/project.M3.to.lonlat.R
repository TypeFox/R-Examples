## ###########################################################
## PURPOSE: Project coordinates from model units to longitude/latitude.
##
## INPUTS:
##   x: x-coordinates of points in model units from project in "file"
##   y: y-coordinates of points in model units from project in "file"
##   file: File name of Models3-formatted file which contains
##         information about the projection (used for "x" and "y").
##   units: Units which x and y are given in.  Options are kilomters ("km")
##     or meters ("m").  The value "deg" does not make sense here, because
##     we would not need to project coordinates from a long/lat reference
##     system to a long/lat reference system.
##   ...:  Other arguments to pass to get.proj.info.M3 function.
##     In this case, the only relevant argument would be the earth
##     radius to use for the projection in "file".
##
## RETURNS: A list containing the elements "coords" and "units".  The
##   element "coords" contains a matrix of coordinates in
##   longitude/latitude.  The element "units" contains the string "deg"
##   to designate that "coords" is in degrees of longitude/latitude.
##
## ASSUMES:
##   1. Availability of R packages ncdf4 and rgdal.
##   2. Projection is lambert conic conformal or polar stereographic. 
##
##
## REVISION HISTORY:
##   Original release: 2011-06-02
## ###########################################################
project.M3.to.lonlat <- function(x, y, file, units, ...){
  
  ## Form projection string describing the projection in the given
  ## Models3-formatted file.
  proj.string <- get.proj.info.M3(file, ...)


  ## Return an error if the Models3-formatted file given is already
  ## referenced by long/lat.  If it is, then we don't need to project
  ## the given coordinates, since they're already on the long/lat system.
  if ( substring(proj.string, first=7, last=13)=="longlat" )
    stop(paste("No need to project, since file ", file,
               " is gridded on long/lat system.", sep=""))

  
  ## ##########################
  ## To avoid errors, the user must specify the units x and y are in.
  ## This must be either "km" or "m".

  if (missing(units))
    stop('User must specify whether units are "km" or "m".')
  else if (units=="km")
    coords.lonlat <- project(1000*cbind(x, y), proj=proj.string, inv=TRUE)
  else if (units=="m")
    coords.lonlat <- project(cbind(x, y), proj=proj.string, inv=TRUE)
  else
    stop(paste(units, " is not a valid option.", sep=""))
  ## ##########################


  ## Name columns appropriately.
  colnames(coords.lonlat) <- c("longitude", "latitude")

  return(list(coords=coords.lonlat, units="deg"))
}
