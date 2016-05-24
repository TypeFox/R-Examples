## ###########################################################
## PURPOSE: Get map lines in the model projection units.
## 
## INPUT:
##   file: File name of Models3-formatted file providing the model
##     projection.
##   database: Geogrpaphical database to use.  Choices are "world",
##     "canusamex", "worldHires", "usa", "state", and "county".
##     Default is "state".
##   units:  Units to be used for the projected coordinates; that is,
##     either "m" (meters) or "km" (kilometers). The option "deg" does
##     not make sense because we would not need to project coordinates
##     from a long/lat reference system to a Models3-formatted file
##     which was also gridded by long/lat.  If unspecified, the default
##     is "km".
##   ...: Other arguments to pass to get.proj.info.M3 function.
##     In this case, the only relevant argument would be the earth
##     radius to use when doing the projections.
##
##
## RETURNS: A list containing two elements "coords" and "units", The
##   element "coords" contains a matrix with the map lines in the
##   projection coordinates.  The element "units" contains the units
##   of the coordinates ("km" or "m").
##
## 
## ASSUMES:
##   Availability of R packages maps, mapdata, ncdf4, and rgdal.  Also
##   uses function get.canusamex.bds in this package.
##
##
## NOTE: The model projection units are of the sort derived in the
##       function proj.lonlat.to.M3().
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-06-02
##   
##   2012-01-11 (JLS): Added additional option of "canmex" for
##     database argument.  Changed argument "region" to argument named
##     "database" for compatibility with the map() function in the "maps"
##     package.
##
##   2012-05-15 (JLS): Changed option "canmex" to "canusmex" and added
##     option "WorldHires".  Both of these use the high-resolution
##     database in package "mapdata" to get better national boundardies
##     (and natural features like coastlines).  The old option "canmex"
##     only used the the relative low-resolution world maps for Canada
##     and Mexico, and was a bit more primitive.
## ###########################################################
get.map.lines.M3.proj <- function(file, database="state", units, ...){

  ## Form projection string describing the projection in the given
  ## Models3 file.
  proj.string <- get.proj.info.M3(file, ...)


  ## Return an error if the Models3 file given is already referenced
  ## by long/lat.  If it is, then we don't need to project the given
  ## coordinates, since they're already on the long/lat system.
  if ( substring(proj.string, first=7, last=13)=="longlat" )
    stop(paste("No need to project, since file ", file,
               " is gridded on long/lat system.", sep=""))


  ## Get the coords of the boundary lines in lat/lon.  If the user
  ## chooses a database like "state", "world", "usa", etc., then we
  ## can get the boundary lines using map() in the "maps" package.  If
  ## the user chooses "canmex", which is not an option for map(), we
  ## call a separate piece of code.
  if (database == "canusamex")
    map.lonlat <- get.canusamex.bds()
  else{
    raw.map.lonlat <- map(database, plot=FALSE, resolution=0)
    map.lonlat <- cbind(raw.map.lonlat$x, raw.map.lonlat$y)
    rm(raw.map.lonlat)
  }


  ## We want to re-project these map boundaries onto the projection that
  ## provided by the specified Models3 file.
  map.CMAQ <- project(map.lonlat, proj=proj.string)
  colnames(map.CMAQ) <- c("x", "y")
  rm(map.lonlat)


  ## ##########################
  ## If the user does not specify desired units, then use "km".  If
  ## user specifies an option other than "km" or "m", give message and
  ## exit function.

  if (missing(units)){
    map.CMAQ <- map.CMAQ/1000
    units <- "km"
  }
  else if (units=="km"){
    map.CMAQ <- map.CMAQ/1000
    units <- "km"
  }
  else if (units=="m")
    units <- "m"
  else
    stop(paste(units, " is not a valid option.", sep=""))
  ## ##########################


  ## Return a list, with the coords in the first position and the
  ## units of those coords in the second.
  x <- list(coords=map.CMAQ, units=units)
  return(x)
}
