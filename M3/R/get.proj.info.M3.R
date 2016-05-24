## ###########################################################
## PURPOSE: Pull information about the projection used from a
##   Models3-formatted file.  Build a string describing the projection
##   which can be used by the R package rgdal.
##
## INPUT:
##   file: File name of Models3-formatted file whose projection we
##     want to use.
##   earth.radius: Assumes a spherical earth, but note that radius may
##     differ in different versions of the Models3 I/O API.  The
##     default is set to the current value (6 370 000 m) in I/O API.
##     An example of another choice is 6 370 997 m, which was the
##     radius used in previous versions of the Models3 I/O API and
##     in previous R packages supplied by Battelle.
##
## RETURNS: String describing model projection, which can be utilized
##   by the rgdal package in R (for projections to and from
##   longitude/latitude, for example).
##
## ASSUMES:
##   1. Your Models3-formatted file uses a Lambert conic conformal or
##      polar stereographic projection, or is gridded in
##      longitude/latitude coordinates.
##   2. Availability of R package ncdf4.
##
## 
## NOTES:
##   1. The R package rgdal depends on the PROJ.4 cartographic
##      projections library (http://trac.osgeo.org/proj/).  The format
##      of the string must therefore be in a style acceptable to
##      PROJ.4.  (See http://www.remotesensing.org/geotiff/proj_list.)
##   2. See information about the meaning of IOAPI projection
##      arguments at
##      http://www.baronams.com/products/ioapi/GRIDS.html.
##   3. Projection info is stored in global attributes of the
##      Models3-formatted file.
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-06-02
## ###########################################################
get.proj.info.M3 <- function(file, earth.radius=6370000){

  ## Open netCDF file which has the projection we want to use..
  nc <- nc_open(file)


  ## ##############################################
  ## Check the validity of earth.radius parameter.  The two most
  ## common entries would probably be 6370000 m (the default) or
  ## 6370997 m (the old Models3 value, as used in the R packages
  ## provided by Battelle).  One could also imagine these values
  ## being entered in kilometers, so we check for that.

  ## Test to see if earth.radius may have been accidentally listed in
  ## kilometers.  If so, it woulld probably be within these bounds.
  if ( (earth.radius > 6360) && (earth.radius < 6380) ){
    warning(paste("Assuming that the given radius, ", earth.radius,
                  ", is in kilomters.  Will use ", earth.radius*1000,
                  "m", sep=""))
    earth.radius <- earth.radius*1000
  }
  else if ((earth.radius < 6360000) || (earth.radius > 6380000))
    warning(paste("Using given radius, ", earth.radius,
                  ", but user should check that this radius is realistic."))
  ## ##############################################
  

  ## ##############################################
  ## FIND OUT THE PROJECTION AND PARAMETERS GOVERNING IT.
  
  ## Find out what projection the grid is on.
  grid.type <- ncatt_get(nc, varid=0, attname="GDTYP")$value

  ## Depending on the type of grid, we extract the information we need
  ## to govern that type of projection.  The projection info is found
  ## in the global attributes.  (To get global attributes, rather than
  ## variable attributes, give 0 as the variable ID.)

   ## Latitude/longitude (if GDTYP==1)
  if (grid.type==1){

    ## Lat/lon projection does not use parameters proj_alpha,
    ## proj_beta, and proj_gamma.
    proj.string <- paste("+proj=longlat", " +a=", earth.radius, " +b=",
                         earth.radius, sep="")
  }


  ## Lambert conic conformal (if GDTYP==2)
  else if (grid.type==2){
    ## Standard parallel 1 is given by P_ALP.
    p.alp <- ncatt_get(nc, varid=0, attname="P_ALP")$value
    ## Standard parallel 2 is given by P_BET.
    p.bet <- ncatt_get(nc, varid=0, attname="P_BET")$value
    ## Central meridian is given by P_GAM.
    p.gam <- ncatt_get(nc, varid=0, attname="P_GAM")$value
    ## Latitude of the center of the Cartesian coordinate system given
    ## by YCENT.
    ycent <- ncatt_get(nc, varid=0, attname="YCENT")$value

    ## Form string based on the projection information.
    proj.string <- paste("+proj=lcc +lat_1=", p.alp, " +lat_2=", p.bet,
                         " +lat_0=", ycent,
                         " +lon_0=", p.gam,
                         " +a=", earth.radius, " +b=", earth.radius, sep="")
  }


  ## Polar stereographic (if GDTYP==6)
  else if (grid.type==6){

    ## P_ALP identifies pole:  north (1) or south (-1)
    p.alp <- ncatt_get(nc, varid=0, attname="P_ALP")$value
    if (p.alp==1.0)
      proj4.lat0 <- 90.0  ##Latitude for North Pole for PROJ.4
    else if (p.alp==-1.0)
      proj4.lat0 <- -90.0  ##Latitude for South Pole for PROJ.4
    else{
      nc <- nc_close(nc)
      stop(paste("For polar stereographic projections (GDTYP=6), P_ALP is ",
                 p.alp, "; it should be either 1 or -1.", sep=""))
    }

    ## P_BET identifies the "secant latitude" (latitude of the true scale).
    p.bet <- ncatt_get(nc, varid=0, attname="P_BET")$value
    ## Central meridian is given by P_GAM.
    p.gam <- ncatt_get(nc, varid=0, attname="P_GAM")$value

    ## Form string based on the projection information.
    proj.string <- paste("+proj=stere +lat_ts=", p.bet,
                         " +lat_0=", proj4.lat0,
                         " +lon_0=", p.gam,
                         " +a=", earth.radius, " +b=", earth.radius, sep="")
  }


  ## Function will be exited with a warning if the grid type is not
  ## one of those listed above.
  else{
    ## Close the Models3 file and exit the function.
    nc <- nc_close(nc)
    stop(paste("Grid type ", grid.type, " cannot be handled by this function.", sep=""))
  }
  ## ##############################################

  
  ## Close the Models3 file.
  nc <- nc_close(nc)
  rm(nc)


  ## Return the string which can be passed to project() and other
  ## functions in R package rgdal.
  return(proj.string)
}
