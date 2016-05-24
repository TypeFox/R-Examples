## ###########################################################
## PURPOSE: For either the rows or the columns, return the coordinates
##   of the centers or the edges of the grid cells.
##
## INPUT:
##   file: File name of Models3-formatted file of interest.
##   dimension: User chooses to get information for either rows
##     ("row") or columns ("column" or "col").
##   position: User chooses whether to get center (default), lower
##     edge, or upper edge coordinates of each row or each column.
##   units: "m" (meters), "km" (kilometers), or "deg" (degrees).  If
##     unspecified, the default is "deg" if the file has a
##     longitude/latitude based grid, and "km" otherwise.
##
## RETURNS: A list containing two elements, "coords" and "units".  If
##   dimension is "row", return as element "coords" a vector
##   containing the y-coordinate of the center, left ("lower"), or
##   right ("upper") edge of each row, depending on the value of
##   argument "position".  If dimension is "column" or "col", return as 
##   element "coords" a vector containing the x-coordinate of the center,
##   left ("lower"), or right ("upper") edge of each row, depending on
##   the value of argument "position".  In both cases, return as element
##   "units" the units of the coordinates (can be "km", "m", or "deg").
##
## ASSUMES:
##   Availability of function get.grid.info.M3().
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-05-20
## ###########################################################
get.coord.for.dimension <- function(file, dimension, position="ctr",
                                    units){

  ## Get info about the grid.
  grid.info <- get.grid.info.M3(file)


  ## Depending on whether we want the center, lower edge, or upper
  ## edge for each cell, we set the offset appropriately.
  if (position=="ctr")
    offset <- 0.5
  else if (position=="lower")
    offset <- 0.0
  else if (position=="upper")
    offset <- 1.0
  else
    stop('Position parameter must be either "lower", "ctr", or "upper".')


  ## Take different actions depending on whether user chooses "row" or
  ## "column" for dimenstions.
  if (dimension=="row")
    coords <- seq(from=grid.info$y.orig+(offset*grid.info$y.cell.width), by=grid.info$y.cell.width, length=grid.info$nrows)
  else if ( (dimension=="column") || (dimension=="col") )
    coords <- seq(from=grid.info$x.orig+(offset*grid.info$x.cell.width), by=grid.info$x.cell.width, length=grid.info$ncols)
  else
    stop('Parameter dimension must be either "row" or "column".')


  ## ##########################
  ## Compare the units returned by get.grid.info.M3() to the
  ## desired units specified by the user.  Convert from m to km, if
  ## necessary.  If the user does not specify desired units, then use
  ## "deg" if the file is based on longitude/latitude and "km"
  ## otherwise.  Warn the user if the units they specified are degrees
  ## ("deg") when the file has meters ("m"), or vice-versa.  If user
  ## specifies an option other than "km", "m", or "deg" for parameter
  ## units, give message and exit function.


  ## If user does not specify units, we look at the units specified in
  ## the object returned by the call to get.grid.info.M3().  If it
  ## it is"deg", we return "deg".  Otherwise we return "km".  
  if (missing(units)){
    if (grid.info$hz.units=="deg")
      units <- "deg"
    else{
      coords <- coords/1000
      units <- "km"
    }
  }

  ## If user specifies "km" we either need to
  ## (1) transform from m, or
  ## (2) warn the user and keep degrees, if grid is in degrees long/lat
  else if (units=="km"){
    if (grid.info$hz.units=="m")  ## Divide by 1000 to go from m to km.
      coords <- coords/1000
    else if (grid.info$hz.units=="deg"){
      warning(paste("Grid specified in file ", file, " is in degrees long/lat; returning degrees.", sep=""))
      units <- "deg"
    }
  }

  ## If the user specifies "m" we just need to make sure the grid is
  ## not in degrees long/lat.  If it is, warn user and keep degrees.
  else if (units=="m"){
    if (grid.info$hz.units=="deg"){
      warning(paste("Grid specified in file ", file, " is in degrees long/lat; returning degrees.", sep=""))
      units <- "deg"
    }
  }

  ## If the user specifies "degrees" we just need to make sure the
  ## grid is specified in degrees long/lat.  If not, we warn user and
  ## return "km".
  else if (units=="deg"){
    if (grid.info$hz.units=="m"){
      warning(paste("Grid specified in file ", file, " is on a projection other than degrees long/lat; returning kilometers.", sep=""))
      coords <- coords/1000
      units <- "km"
    }
  }

  ## If the user specifies something other than "m", "km", or "deg",
  ## we stop and print an error message.
  else
    stop(paste(units, " is not a valid option.for 'units'.", sep=""))
  ## ##########################


  ## Return a list, with the coords in the first position and the
  ## units of those coords in the second.
  x <- list(coords=coords, units=units)
  return(x)
}
