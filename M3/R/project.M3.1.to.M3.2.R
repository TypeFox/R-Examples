## ###########################################################
## PURPOSE: Project coordinates based on projection in first
##   Models3-formatted file to the projection given in second
##   Models3-formatted file.
##
## INPUTS:
##   x: x-coordinates in model units from projection in from.file
##   y: y-coordinates in model units from projection in from.file
##   from.file: Name of Models3-formatted file with the model
##     projection underlying "x" and "y"
##   to.file: Name of Models3-formatted file with the model
##     projection to which you want "x" and "y" to be projected.
##   units: Units of "x" and "y".  The coordinates returned will also
##     be in these units.
##   ...: Other arguments to pass to get.proj.info.M3 function.
##     In this case, the only relevant argument would be the earth
##     radius to use when doing the projections.
##
## RETURNS: A list containing two elements "coords" and "units".  The
##   element "coords" contains a matrix of coordinates using projection
##   in to.file.  The element "units" contains the units of the
##   coordinates, which are the same as those specified for input "x"
##   and "y".
##
## ASSUMES:
##   1. Availability of R packages ncdf4 and rgdal.
##   2. Projections are lambert conic conformal or polar stereographic.
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-06-02
## ###########################################################
project.M3.1.to.M3.2 <- function(x, y, from.file, to.file, units, ...){

  ## Form projection string describing the projection in the given
  ## Models3-formatted file.
  from.proj.string <- get.proj.info.M3(from.file, ...)
  to.proj.string <- get.proj.info.M3(to.file, ...)



  ## ##########################
  ## If both files have the same projection, then we don't need to do any
  ## projections.  Exit, and tell user this.
  
  ## Return an error if both files are already on the same projection.
  if ( identical(from.proj.string, to.proj.string) )
    stop(paste("No need to project, since file ", from.file,
               " and file ", to.file, " use the same projection.", sep=""))
  ## ##########################



  ## ########################## 
  ## Take different course if one of the projections is long/lat.
  
  ## (1) If "from" projection is long/lat, then we need function
  ## project.lonlat.to.M3.  Warn user that this funciton is
  ## being called.
  if ( substring(from.proj.string, first=7, last=13)=="longlat" ){
    warning(paste('The specified "from" projection in file,', from.file, ' is longitude/latitude.  Using function project.lonlat.to.M3...', sep=""))
    return(project.lonlat.to.M3(longitude=x, latitude=y, file=to.file,
                                units=units, ...))
  }

  
  ## (2) If "to" projection is long/lat, then we need function
  ## project.M3.to.lonlat(). Warn user that this funciton is
  ## being called.
  if ( substring(to.proj.string, first=7, last=13)=="longlat" ){
    warning(paste('The specified "to" projection in file,', to.file, ' is longitude/latitude.  Using function project.M3.to.lonlat.', sep=""))
    return(project.M3.to.lonlat(x=x, y=y, file=from.file, units=units, ...))
  }

  
  ## (3) If neither project is long/lat, proceed with the remainder of
  ## this function.

  ## ##########
  ## Put the given coordinates into "Spatial Points" form for use
  ## with spTransform function.  Make sure we adjust properly for
  ## units (must be either "km" or "m".)
  
  if (missing(units))
    stop('User must specify units of "x" and "y".')
  else if (units=="km")
    from.coords <- data.frame(x=1000*x, y=1000*y)
  else if (units=="m")
    from.coords <- data.frame(x=x, y=y)
  else
    stop(paste(units, " is not a valid option for 'units'.", sep=""))
  ## ##########


  ## ##########
  ## Make into a SpatialPoints object.
  coordinates(from.coords) <- c("x", "y")
  proj4string(from.coords) <- CRS(from.proj.string)
  
  ## Project locations from CMAQ units to longitude and latitude.
  to.coords <- spTransform(from.coords, CRS(to.proj.string))
  
  ## Return the coordinates, adjusting units appropriately 
  if (units=="km")
    return(list(coords=coordinates(to.coords)/1000, units=units))
  else
    return(list(coords=coordinates(to.coords), units=units))
  ## ##########
  
  ## ##########################
}
