## ###########################################################
## PURPOSE: Obtain outlines of Canada, USA, and Mexico, including
## state lines.
##
##
## INPUT: No input arguments.
##
## 
## RETURNS: A list containing two elements "coords" and "units", The
##   element "coords" contains a matrix with the map lines in the
##   projection coordinates.  The element "units" contains the units
##   of the coordinates ("km" or "m").
##
##
## ASSUMES:
##   Availability of R packages maps and mapdata.
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2012-05-15
## ###########################################################
get.canusamex.bds <- function(){

  ## Get maps of Canada, USA, Mexico in high-resolution from mapdata
  ## package.
  canusamex.natl <- map("worldHires", regions=c("Canada", "USA", "Mexico"),
                   exact=FALSE, resolution=0, plot=FALSE)
  canusamex.natl.lonlat <- cbind(canusamex.natl$x, canusamex.natl$y)
  rm(canusamex.natl)


  ## Get map of state borders, without including the outer boundaries.
  ## These outer boundaries are provided by the the high-resolution
  ## national maps created in the previous step.
  state <- map("state", exact=F, boundary=FALSE, resolution=0, plot=FALSE)
  state.lonlat <- cbind(state$x, state$y)
  rm(state)


  ## Put it national and state boundaries together.
  canusamex <- rbind(canusamex.natl.lonlat, matrix(NA, ncol=2), state.lonlat)
  rm(canusamex.natl.lonlat, state.lonlat)


  ## Return the boundaries to the calling program/function.
  return(canusamex)
}
