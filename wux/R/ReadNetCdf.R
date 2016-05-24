
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2014-09-29 11:36:55 +0200 (Mon, 29 Sep 2014) $
# $Rev: 286 $
# ----------------------------------------------------------------

## ----------------------------------------------------------------
## Input routines and processes handling input routines live here:
##
## - Process: get aggregated spatio-temporal data
## - Reading NetCDF data arrays, and grid arrays
## - appertaining NetCDF helper funtions
## - Filename handling
## ----------------------------------------------------------------

ReadNetCdf <- function(filenames, parameter.name,
                       time.vectors, adjust.values = TRUE,
                       offset=FALSE, count=FALSE,
                       indices.exclude=FALSE, ...) {
  ## Reads in all NetCDF files specified by "filenames" and concatenates them.
  ##
  ## Args:
  ##   filenames: Character vector with filenames to be read in.
  ##   parameter.name: Variable name of parameter in NetCDF file.
  ##   time.vectors: Time structure returned from ReadNetCdfTime.
  ##   adjust.values: Should the data values units be transformed?
  ##   offset, count: NetCDF file offset and count.
  ##   indices.exclude: Matrix indicating what values will be excluded from
  ##                    data array.
  ##   '...': Currently no additional meaning.
  ##
  ## Returns:
  ##   3-Dimensional array with data values,
  ##   optionally not whole file, but only subregion given by 'offset',
  ##   'count' and 'index.exclude'.
  ##
  ## History:
  ##   2010-10-25 | Original code (thm)
  ##   2011-02-04 | missing values and units extracted now (geh, thm)
  ##   2011-05-26 | reads in data with dim = 1 as well (thm)
  ##   2011-06-01 | reads in data with dim = 2 as well (thm)
  ##   2011-12-01 | some small renaming (thm)
  ##   2014-09-29 | _FillValue attribute is also affected by scaling and offset... (thm)
  ##   2016-01-13 | change to "ncdf4" library

  ## get parameter longname (CF-convention)
  par.longname <- names(parameter.name)

  ## if more datanames have been passed on
  for (infile in filenames) {
    ## get filename only (without whole path)
    file.name <- tail(unlist(strsplit(infile, .Platform$file.sep)), n=1)
    cat("\n      opening file \"", file.name, "\"\n", sep="")

    ## open nc...
    ## nc <- open.ncdf(infile, readunlim=TRUE)
    nc <- nc_open(infile, readunlim=TRUE)

    ## Does our specified variable exist?
    if ( !is.list(parameter.name) ) {
      if ( is.null(nc$var[[parameter.name]]) )
        stop("    ERROR: SHORT PARAMETER NAME NOT LOCATED AS VARIABLE ",
             "NAME IN NC-FILE")
    } else {
      nc.varcheck.names <- names(nc$var[parameter.name[[1]]])
      var.pos <- which(!is.na(nc.varcheck.names))
      if ( length(var.pos) == 1 ) {
        parameter.name <- nc.varcheck.names[var.pos]
      } else if ( length(var.pos) ==  0 ) {
        stop("    ERROR: SHORT PARAMETER NAME NOT LOCATED AS VARIABLE ",
             "NAME IN NC-FILE")
      } else if ( length(var.pos) >=  1 ) {
        stop("    ERROR: MULTIPLE SHORT PARAMETER NAMES ARE LOCATED ",
             "AS VARIABLE NAME IN NC-FILE")
      }
    }

    ## read dimensions of netcdf file
    nc.dims <- GetNetcdfDims(nc)

    ## read dimensions of parameter
    parameter.dims <- GetVariableDims(nc, parameter.name)

    ## setting offset default 1
    nc.start <- rep(1, times=length(parameter.dims))
    names(nc.start) <- names(parameter.dims)
    ## paste lon lat offset of given 'offset' argument
    if (length(offset) >= 2) {
      nc.start[c("lon", "lat")] <- offset
      nc.start[c("time")] <- time.vectors[[infile]][["offset"]]
    }
    ##read in whole dimension range for count
    nc.count <- parameter.dims
    ## paste lon lat counts of passed 'count' argument
    if (length(count) >= 2) {
      nc.count[c("lon", "lat")] <- count
      nc.count[c("time")] <-  time.vectors[[infile]][["count"]]
    }

    ## addoffset and scale factor for packed data. if data
    ## are not packed functions return add.offset = 0 and
    ## scale.factor = 1
    ## see http://www.unidata.ucar.edu/software/netcdf/docs/BestPractices.html
    add.offset <- GetAddOffset(nc, parameter.name)
    scale.factor <- GetScaleFactor(nc, parameter.name)

    ## write data to memory
    time.start <- as.character(time.vectors[[infile]][["intervall.begin"]])
    time.end <- as.character(time.vectors[[infile]][["intervall.end"]])

    cat("      count:", names(nc.count), "\n")
    cat("             ", nc.count, "\n")

    cat("      reading parameter data \"", parameter.name, "\" from ",
        time.start, " to ", time.end, "\n", sep = "")
    ## old ncdf package...
    ## nc.values <- get.var.ncdf(nc, parameter.name,
    ##                           start=nc.start, count=nc.count, verbose = FALSE)
    nc.values <- ncvar_get(nc, parameter.name,
                              start=nc.start, count=nc.count, verbose = FALSE)

    ## force nc.values to be a 3-dimensional array (as nc.count)
    n.dims <- length(dim(nc.values))
    if (n.dims == 1)
     dim(nc.values) <- nc.count[c("lon", "lat", "time")]
    else if (n.dims == 2)
      dim(nc.values) <- nc.count[c("lon", "lat", "time")]
    else if (n.dims > 3)
      stop("HOW DID YOU GET HERE? this has to be implemented...")
    ## maybe generalize dim extension by something like this:
    ## nc2 <- (nc.values)
    ## dim(nc2) <- c(lat = 1, dim(nc2))
    ## expand.dim <- !names(nc.count[c("lon","lat")]) %in% names(dim(nc2))
    ## dim(nc2) <- c(nc.count[c("lon","lat")][expand.dim], dim(nc2))

    ## set missing values to NA
    ## missing.list <- att.get.ncdf(nc, parameter.name, "_FillValue")
    missing.list <- ncatt_get(nc, parameter.name, "_FillValue")
    if (missing.list$hasatt) {
      missing.value <- missing.list$value * scale.factor + add.offset
      nc.values[nc.values == missing.value] <- NA
    }

    ## check out units
    ## units <- att.get.ncdf(nc, parameter.name, "units")
    units <- ncatt_get(nc, parameter.name, "units")
    if (units$hasatt){
      units.value <- units$value
    } else {
      stop("PARAMETER UNITS NOT DECLARED! PASS THEM THROUGH USER.INPUT SOMEHOW")
    }

    ## close.ncdf(nc)
    nc_close(nc)

    ## retreiving dimension of output array
    nc.array.dim <- attr(nc.values, "dim")
    last.dim <- length(nc.array.dim)
    ## For every timestep (i.e. over dim 3) set values outside subregion to NA.
    ## The array function with dim statement is executed to keep the array as it
    ## has been read in.
    cat("      clipping subregions\n")
    nc.values <- array(apply(nc.values, 3,
                             function(x){is.na(x) <- indices.exclude; x}),
                       dim = nc.array.dim)

    ## append array to previous array
    if (length(ls(pattern="nc.values.concatenated")) == 0) {
      ## no data appended yet
      nc.values.concatenated <- nc.values
    }  else  {
      cat("\n      concatenating fields\n")
      ## appending
      nc.values.concatenated <- abind(nc.values.concatenated, nc.values)
      dimnames(nc.values.concatenated) <- NULL
    }
    
  }
  ## adjust values, i.e. bring to common scale Celsius, mm...
  if (adjust.values == TRUE)
    nc.values.concatenated <- AdjustValues(nc.values.concatenated,
                                         par.longname, units.value)

  return(nc.values.concatenated)
}



AdjustValues <- function(x, parameter.name, units){
  ## Transforms units. Kelvin to Celsius, kg/(m2*s) to mm/day.
  ##
  ## Args:
  ##   x: data to be transformed.
  ##   parameter.name: Name of parameter corresponding to CF-conventions.
  ##   units: Character. Units of given parameter.
  ##
  ## Returns:
  ##   Transformed data array with new units.
  ##
  ## History:
  ##   2010-10-27 | Original code (thm)
  ##   2011-02-04 | improved units handling (geh, thm)
  ##   2012-01-09 | generalized air temperature case (for min and max) (thm)

  ## get rid of space
  units <- stringr::str_trim(units)

  if (parameter.name == "precipitation_amount") {
    kg.units <- c("kgm-2s-1", "kg m-2 s-1")
    mm.units <- c("mm", "mm/day", "kg m-2", "mm/d", "kg m-2 d-1", "kg/m^2")
    m.units <- c("m")
    if (units %in% c(kg.units, mm.units, m.units)) {
      if (units %in% mm.units)
        x <- x
      if (units %in% kg.units)
        x <- x * 86400
      if (units %in% m.units)
        x <- x * 1000
    } else {
      stop("UNKNOWN PARAMETER UNIT. PLEASE ADD IT TO AdjustValues()")
    }
  }

  if (length(grep("air_temperature", parameter.name)) > 0) {
    if (units %in% c("K", "C", "degree C", "Celsius")) {
      if (units == "K")
        x <- x - 273.15
      if (units %in% c("C", "C", "degree C", "Celsius"))
        x <- x
    } else {
      stop("UNKNOWN PARAMETER UNIT. PLEASE ADD IT TO AdjustValues()")
    }
  }

  return(x)
}


GetNetcdfDims <- function(nc) {
  ## Extracts dimension information from nc object.
  ##
  ## Args:
  ##   nc: Class 'ncdf' or 'RasterLayer' (currently not in use) object.
  ##
  ## Returns:
  ##   Vector containing dimensions from nc-object (NetCDF file), named
  ##   'lon', 'lat', or 'time' (e.g. instead of 'x', 'y', 'TIME').
  ##
  ## History:
  ##   2010-10-29 | Original code (thm)
  
  ## if nc of 'netcdf' package
  if (class(nc) %in% c('ncdf', 'ncdf4')) {
    nc.dims <- rep(0, times=nc$ndims)
    ## give nc.dims the netCDF dimension names and dimension lengths
    for (ii in seq(along = nc$dim)) {
      ## get NetCDF dimension name 
      names(nc.dims)[ii] <- nc$dim[[ii]]$name
      ## get dimension value
      nc.dims[ii] <- nc$dim[[ii]]$len
    }
  }
  ## if nc of 'raster' package class RasterLayer
  if (nc$ndims == 3) {
    if (class(nc) == 'RasterLayer') {
      nc.dims[1] <- nc@nrows
      nc.dims[2] <- nc@ncols
      ## only one timestep possible in RasterLayer object
      nc.dims[3] <- 1
      names(nc.dims) <- c("lon", "lat", "time")        
    }
  } else {}

  ## name dimanesion names properly (i.e. to 'lon, 'lat' and so on)
  nc.dims <- RenameDimensionNames(nc.dims)
  
  return(nc.dims)
}

GetVariableDims <- function(nc, parameter) {
  ## Extracts dimensions for current 'parameter' variable from nc object.
  ##
  ## Args:
  ##   nc: 'ncdf' class object.
  ##   parameter: Name of variable in nc object (i.e. NetCDF file)
  ##
  ## Returns:
  ##   Named vector with variable dimensions from NetCDF file. Dimension
  ##   names will be renamed to WUX style ('lon', 'lat' and 'time').
  ## 
  ## History:
  ##   2010-10-29 | Original code (thm)

  ## if nc of 'netcdf' package
  if (class(nc) %in% c('ncdf', 'ncdf4')) {
    parameter.var <- nc$var[[parameter]]
    if (is.null(parameter.var))
      stop("SHORTNAME OF PARAMETER IN FILE NOT KNOWN")
    parameter.ndims <- parameter.var$ndims
       
    ## retrieve variable dimensions
    parameter.dims <- parameter.var$size  # pasting vardims
    ## retrieve variable dimension names
    for (ii in seq(parameter.ndims))
      names(parameter.dims)[ii] <- parameter.var$dim[[ii]]$name
 
  } else {
    stop("UNKNWON NC CLASS")
  }
  ## rename to dimension names WUX is working with
  parameter.dims <- RenameDimensionNames(parameter.dims)
  return(parameter.dims)
}

RenameDimensionNames <- function(x) {
  ## Replaces 'names' attribute of vector x by dimension names used
  ##  in WUX ('lon, 'lat' and 'time').
  ## e.g.
  ##   bnds   rlon   rlat height   time 
  ##     2    170    190      1   3600
  ##   to
  ##   bnds   lon   lat height   time 
  ##      2    170   190      1   3600
  ##
  ## Args:
  ##   x: Named dimension vector.
  ##
  ## Returns:
  ##   Same vector as input x, but with replaced 'names' attribute.
  ##
  ## History:
  ##   2010-10-29 | Original code (thm)

  lat.names <- c("y", "rlat", "lat", "latitude", "nlat", "ny", "grid_latitude")
  lon.names <- c("x", "rlon","lon","longitude","nlon", "nx", "grid_longitude")
  time.names <- c("time", "TIME")
  ##  bnds.names <- c("bnds", "time_bnds")
  ## 'names' attribute of known.names are the names to be replaced by
  known.names <- list("lat" = lat.names, "lon" = lon.names, "time" = time.names)

  ## Replace
  for (ii in names(known.names)) {
    ## look for match in known.means values vector
    m <- match(known.names[[ii]], names(x))
    ## getting rid of NAs
    ## in case all is NA (menaing unknown names):
    ## if(all(is.na(m)))
    ##   stop("UNKNOWN DIMESION NAME FOR ", ii, " IN NETCDF FILE. CHECK IN RenameDimensionNames FOR NAMING OR GIVE YOUR NETCDF FILE PROPER DIMENSION NAMES ;) ")
    ## else
      m <- m[!is.na(m)]  

    ## change name in netcdf vector with name of known.names vector
    names(x)[m] <- ii
  }
  return(x)    
}

GetAddOffset <- function(nc, parameter) {
  ## Gets the Add Offset value for unpacking nc data
  ## We still need Scale Factor (see "GetScaleFactor").
  ## Actually there is no need for this function, as the 'ncdf'
  ## input routine is smart enough to unpack automatically.
  ## 
  ## Returns:
  ##   addOffset for packed data (hasAddOffset) or 0 unpacked.
  ## 
  ## History:
  ##   2010-10-29 | Original code (thm)

  ## if nc of 'netcdf' package
  if (class(nc) %in% c('ncdf', 'ncdf4')) {
    ## data packed?
    has.add.offset <- nc$var[[parameter]]$hasAddOffset
    if (is.null(has.add.offset)) {
      ## adding with 0.. nothing happens :)
      add.offset <- 0
    } else {
      if (has.add.offset == TRUE) {
        add.offset <- nc$var[[parameter]]$addOffset
        cat("add.offset = ",add.offset,"\n")
      } else {
        add.offset <- 0
      }
    }      
  }
  return(add.offset)
}

GetScaleFactor <- function(nc, parameter) {
  ## Gets the Scale Factor value for unpacking nc data
  ## Actually there is no need for this function, as the 'ncdf'
  ## input routine is smart enough to unpack automatically.
  ## 
  ## Returns:
  ##    Scale factor for packed data (hasScaleFactor) or 1 for unpacked.
  ## 
  ## History:
  ##   2010-10-29 | Original code (thm)

  ## if nc of 'netcdf' package
  if (class(nc) %in% c('ncdf', 'ncdf4')) {
    ## data packed?
    has.scale.factor <- nc$var[[parameter]]$hasScaleFact
    if (is.null(has.scale.factor)) {
      ## multiplying with 1.. nothing happens
      scale.factor <- 1
    } else {
      if (has.scale.factor == TRUE) {
        scale.factor <- nc$var[[parameter]]$scaleFact
        cat("scale.factor = ",scale.factor,"\n")
      } else {
        scale.factor <- 1
      }
    }
  }      
  return(scale.factor)
}
