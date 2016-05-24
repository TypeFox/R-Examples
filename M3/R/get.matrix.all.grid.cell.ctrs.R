## ###########################################################
## PURPOSE: Find the locations of the grid cell centers in grid
##   units.
##
## INPUT:
##   file: File name of Models3-formatted file of interest.
##   units: "m" (meters), "km" (kilometers), or "deg" (degrees).  If
##     unspecified, the default is "deg" if the file has a
##     longitude/latitude based grid, and "km" otherwise.
##
## RETURNS: A list containing two elements "coords" and "units".  The
##   element "coords" contains a matrix with number of rows equal to the
##   number of grid cells and two columns.  The first column contains the
##   x-coordinate of the grid cell centers; the second column contains
##   the y-coordinate of the grid cell centers.  The points are
##   listed in order such that the x-coordinates are changing faster
##   than the y-coordinates.  The element "units" contains the units of
##   the coordinates (can be "km", "m", or "deg").
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-05-19
## ###########################################################
get.matrix.all.grid.cell.ctrs <- function(file, units){

  ## If user specifies the desired units, pass those on to the
  ## function get.coord.for.dimension().
  if (!missing(units)){
    x.coord <- get.coord.for.dimension(file, dimension="column",
                                       position="ctr", units=units)
    y.coord <- get.coord.for.dimension(file, dimension="row",
                                       position="ctr", units=units)
  }
  ## If user doesn't specify units, don't try to pass them.  The
  ## function get.coord.for.dimension() will return "deg" if the file
  ## is long/lat, and "km" otherwise.
  else{
    x.coord <- get.coord.for.dimension(file, dimension="column",position="ctr")
    y.coord <- get.coord.for.dimension(file, dimension="row", position="ctr")
  }
    
  ## The units for the x and y-coordinates should be the same.  If
  ## not, something very strange has happened.
  if (!identical(x.coord$units, y.coord$units))
    stop(paste("Error: Units for x-coordinates and y-coordinates differ.  For x, units are ", x.coord$units, "; for y, ", y.coord$units, ".", sep=""))


  ## Return a list with the first element a matrix with all combinations
  ## of these x- and y-coordinates and the units as the second element.
  all.ctrs <- as.matrix(expand.grid(x.coord$coords, y.coord$coords))
  colnames(all.ctrs) <- c("x", "y")
  
  return(list(coords=all.ctrs, units=x.coord$units))
}
