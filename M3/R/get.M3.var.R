## ###########################################################
## PURPOSE: Read in variable values from Models3-formatted files.
##
## INPUT:
##   file: name of Models3-formatted file to be read
##   var: name (character string) or number (positive integer) of
##     variable whose values are to be read
##   lcol, ucol: Lower and upper column bounds (positive integers) to
##     be read.  The default is to read all columns.
##   lrow, urow: Lower and upper row bounds (positive integers) to be
##     read.  The default is to read all rows.
##   llay, ulay: Lower and upper layer bounds (positive integers) to
##     be read.  The default is to read the first layer only.
##   ldatetime, udatetime: Starting and ending date-times (either Date
##     or POSIX class) in GMT.  The default is to read all date-times.
##     If the file is time-independent, the one available time step
##     will be read and user input for ldatetime and udatetime will be
##     disregarded.
##   hz.units: Units associated with grid cell horizontal dimensions.
##     Default is degrees ("deg") if the data is indexed according to
##     longitude/latitude and kilometers ("km") otherwise.  If the
##     file is not indexed according to longitude/latitude, the user
##     could specify meters ("m").
##
##
## RETURNS: List with several elements.
##   Element "data" holds the actual variable values in a 4D (or 3D,
##     in the case of time-independent files) array.  Dimensions are
##     columns, rows, layers, date-time steps.
##   Element "data.units" holds the units associated with "data".
##   Elements "x.cell.ctr" and "y.cell.ctr" give the coordinates of the
##     centers of the grid cells in units given by "hz.units".
##   Element "hz.units" gives the units associated with the
##     horizontal dimensions of the grid cells (km, m, or deg).
##   Element "rows" gives the row numbers extracted.
##   Element "cols" gives the column numbers extracted.
##   Element "layers" gives the layer numbers extracted.
##   Element "datetime" gives the date-time steps associated with the
##     variable, if the file is not time-independent.  For
##     time-independent files, element datetime is set to NULL.
##
## ASSUMES: The Models3-formatted file is either time-independent or
##   time-stepped.  It cannot be of type circular-buffer.
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-06-02
## ###########################################################
get.M3.var <- function(file, var, lcol, ucol, lrow, urow,
                            llay, ulay, ldatetime, udatetime,
                            hz.units){

  ## Open netCDF file which has the projection we want to use..
  nc <- nc_open(file)


  ## ##########################
  ## MAKE SURE THE VARIABLE SPECIFIED IN PARAMETER var IS VALID.

  ## Get list of variable names.  Check that the variable name
  ## provided is on the list.  If a variable number is provided, make
  ## sure that there is a variable with that number.
  all.varnames <-  names(nc$var)
  

  ## Check to make sure the user specified a variable to read.  If
  ## not, close file, print error message, and exit.
  if (missing(var)){
    nc <- nc_close(nc)
    stop( paste("Specify the name or number of the variable to be read.  Variable names are: ", paste(all.varnames, collapse=", "), sep="") )
  }


  ## If var is a character string, its name must be on the list of
  ## variable names.
  if ( is.character(var) && !(var %in% all.varnames) ){
    nc <- nc_close(nc)
    stop( paste("File ", file, " does not contain variable named ", var, sep="") )
  }


  ## If var is numeric, it must be an integer between 1 and the number
  ## of variables.
  if (is.numeric(var)){

    ## Specified variable number must be an integer.
    if ( !(trunc(var) == var) ){
      nc <- nc_close(nc)
      stop( paste("Parameter var must be a whole number between 1-", length(all.varnames), sep="") )
    }

    ## Specified variable number cannot be less than 1 or more than
    ## the number of variables.    
    if ( (var < 1) || (var > length(all.varnames)) ){
      nc <- nc_close(nc)
      stop( paste("File ", file, " contains variables numbered 1-", length(all.varnames), sep="") )
    }
  }


  ## Parameter var must be either numeric or a character string.
  if ( !is.numeric(var) && !is.character(var) ){
    nc <- nc_close(nc)
    stop( "Parameter var must give either the name or number of the variable to be read." )
  }
  ## ##########################


  ## ##########################
  ## MAKE SURE THE INPUT FOR WHICH ROWS, COLUMNS, AND LAYERS TO READ
  ## MAKES SENSE.

  ## Find out the dimensions of the chosen variable.  I assume that
  ## these are listed in terms of number of columns, number of rows,
  ## number of layers, number of date-time steps.
  dimens <- nc[["var"]][[var]][["size"]]
  if (length(dimens) < 4){
    nc <- nc_close(nc)
    stop( "There are less than 4 dimensions in this file.")
  }


  ## If lower/upper column limits are missing, then make them the
  ## minimum/maximum available in the file.
  if (missing(lcol))
    lcol <- 1
  if (missing(ucol))
    ucol <- dimens[1]
  ## Check to make sure that the upper column limit is greater than or
  ## equal to the lower column limit.
  if (ucol < lcol){
    nc <- nc_close(nc)
    stop(paste("Upper column limit, ", ucol, ", is less than lower column limit, ", lcol, sep=""))
  }

  
  ## If lower/upper row limits are missing, then make them the
  ## minimum/maximum available in the file.
  if (missing(lrow))
    lrow <- 1
  if (missing(urow))
    urow <- dimens[2]
  ## Check to make sure that the upper row limit is greater than or
  ## equal to the lower row limit.
  if (urow < lrow){
    nc <- nc_close(nc)
    stop(paste("Upper row limit, ", urow, ", is less than lower column limit, ", lrow, sep=""))
  }

  
  ## If lower/upper row limits are missing, then make them the
  ## minimum/maximum available in the file.
  if (missing(llay))
    llay <- 1
  if (missing(ulay))
    ulay <- dimens[3]
  ## Check to make sure that the upper layer limit is greater than or
  ## equal to the lower layer limit.
  if (ulay < llay){
    nc <- nc_close(nc)
    stop( paste("Upper layer limit, ", ulay, ", is less than lower layer limit, ", llay, sep="") )
  }


  ## Check to make sure row, column, and layer limits are positive numbers.
  if ( (lcol <= 0) && (ucol <= 0) && (lrow <= 0) && (urow <= 0) && (llay <= 0) && (ulay <= 0)){
    nc <- nc_close(nc)
    stop("Upper and lower row, column, and layer boundaries must be positive.")
  }


  ## Find the proper subset of columns.
  col.seq <- 1:dimens[1]
  which.col <- which( (lcol <= col.seq) & (col.seq <= ucol) )
  start.col <- min(which.col)
  count.col <- max(which.col) - min(which.col) + 1
  
  ## Find the proper subset of rows.
  row.seq <- 1:dimens[2]
  which.row <- which( (lrow <= row.seq) & (row.seq <= urow) )
  start.row <- min(which.row)
  count.row <- max(which.row) - min(which.row) + 1
  
  ## Find the proper subset of layers.
  lay.seq <- 1:dimens[3]
  which.lay <- which( (llay <= lay.seq) & (lay.seq <= ulay) )
  start.lay <- min(which.lay)
  count.lay <- max(which.lay) - min(which.lay) + 1
  ## ##########################


  ## ##########################
  ## SUBSET THE DATE-TIMES (parameters ldatetime, udatetime).

  ## Check whether the file is time independent.  It is time
  ## independent if the time step increment is zero.
  tstep.incr <- ncatt_get(nc, varid=0, attname="TSTEP")$value


  ## If the time step is not 0, ensure we get the correct range of
  ## time steps.
  if (tstep.incr != 0){

    ## Form a sequence of all the datetimes included in the Models3 file.
    datetime.seq <- get.datetime.seq(file)


    ## If ldatetime is missing, assign it the earliest date-time; if
    ## udatetime is missing assign it the latest date-time.
    if (missing(ldatetime))
      ldatetime <- min(datetime.seq)
    if (missing(udatetime))
      udatetime <- max(datetime.seq)
  

    ## Check to see if the date-time limits are in Date format.  If so,
    ## make them into a POSIX format date.  For the lower limit, this
    ## would mean a time stamp at midnight (beginning of the given
    ## date).  For the upper limit, this would mean a time stamp at 23:59:59
    ## (last part of the given date).
    if ("Date" %in% class(ldatetime))
      ldatetime <- combine.date.and.time(date=ldatetime, time="00:00:00")
    if ("Date" %in% class(udatetime))
      udatetime <- combine.date.and.time(date=udatetime, time="23:59:59")


    ## Check to make sure lower bound on datetime is same as or earlier
    ## than the upper bound.
    if (udatetime < ldatetime){
      nc <- nc_close(nc)
      stop(paste("Upper date-time bound, ", udatetime, ", is before lower date-time bound, ", ldatetime, sep=""))
    }


    ## Find the dates in the sequence which fall in the specified range.
    which.datetime <- which( (ldatetime <= datetime.seq) & (datetime.seq <= udatetime) )
    start.datetime <- min(which.datetime)
    count.datetime <- max(which.datetime) - min(which.datetime) + 1
  }


  ## For a time indep. file, can only read the one time step available.
  else{
    message("Time independent file - reading only time step available.")
    start.datetime <- 1
    count.datetime <- 1
  }
  ## ##########################


  ## ##########################
  ## ACTUALLY EXTRACT THE DATA FOR THIS VARIABLE (parameter var).

  extracted.data <- ncvar_get( nc, varid=var, start=c(start.col, start.row, start.lay, start.datetime), count=c(count.col, count.row, count.lay, count.datetime) )

  ## If the time step is not 0, then there are presumably meaningful
  ## time steps in this file.  If the time step is 0, then the file is
  ## time-independent, and we allow the array to be a 3D array.
  
  ## Force this matrix/array into a 4D array, so that all dimensions
  ## are represented, even if some dimensions are of length 1.  Then,
  ## we pass the units for each dimension.
  if (tstep.incr != 0)
    dim(extracted.data) <- c(count.col, count.row, count.lay,
                           count.datetime)
  else
    dim(extracted.data) <- c(count.col, count.row, count.lay)


  ## If units for this data given, store them.  Otherwise, mark the
  ## data units as missing.
  info.data.units <-  ncatt_get(nc, varid=var, attname="units")
  if (info.data.units$hasatt)
    data.units <- info.data.units$value
  else
    data.units <- NA
  rm(info.data.units)
  ## ##########################

  
  ## ##########################
  ## FIND THE UNITS ASSOCIATED WITH X- AND Y- COORDINATES (IN MODEL
  ## UNITS) OF THE CENTER OF THE GRID CELLS.  IF USER HAS SPECIFIED
  ## THE DESIRED UNITS, WE PASS THAT ALONG TO
  ## get.coord.for.dimension(); OTHERWISE, WE TAKE DEFAULT UNITS
  ## RETURNED BY get.coord.for.dimension().  WE ALREADY HAVE THE
  ## DATETIME SEQUENCE FROM THE CALCULATIONS ABOVE.
  if (missing(hz.units)){
    all.x.coord <- get.coord.for.dimension(file, dimension="col",
                                       position="ctr")
    all.y.coord <- get.coord.for.dimension(file, dimension="row",
                                       position="ctr")
  }
  else{
    all.x.coord <- get.coord.for.dimension(file, dimension="col",
                                       position="ctr", units=hz.units)
    all.y.coord <- get.coord.for.dimension(file, dimension="row",
                                       position="ctr", units=hz.units)
  }

  ## The units for the x and y-coordinates should be the same.  If
  ## not, something very strange has happened.
  if (!identical(all.x.coord$units, all.y.coord$units))
    stop(paste("Error: Units for x-coordinates and y-coordinates differ.  For x, units are ", all.x.coord$units, "; for y, ", all.y.coord$units, ".", sep=""))
  else
    hz.units <- all.x.coord$units
  

  ## If the user has specified only a certain subset of rows and
  ## columns to be returned, we subset the coordinates appropriately.
  x.coord <- all.x.coord$coords[which.col]
  y.coord <- all.y.coord$coords[which.row]
  rm(all.x.coord, all.y.coord)
  ## ##########################


  ## ##########################
  ## PUT DATA AND UNITS TOGETHER IN A LIST TO RETURN TO USER.

  ## If not a time-independent file, then include date-time steps.
  if (tstep.incr != 0)
    extracted.list <- list(data=extracted.data, data.units=data.units,
                           x.cell.ctr=x.coord, y.cell.ctr=y.coord,
                           hz.units=hz.units,
                           rows=row.seq[which.row],
                           cols=col.seq[which.col],
                           layers=lay.seq[which.lay],
                           datetime=datetime.seq[which.datetime])
  else
    extracted.list <- list(data=extracted.data, data.units=data.units,
                           x.cell.ctr=x.coord, y.cell.ctr=y.coord,
                           hz.units=hz.units,
                           rows=row.seq[which.row],
                           cols=col.seq[which.col],
                           layers=lay.seq[which.lay],
                           datetime=NULL)
  ## ##########################

  
  ## Close netCDF file.
  nc <- nc_close(nc)
  
  ## Return list of extracted information about variable var.
  return(extracted.list)
}
