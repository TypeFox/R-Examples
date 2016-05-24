## ###########################################################
## PURPOSE: Subset the list resulting from a get.M3.var()
##   function call.
##
## INPUT:
##   var.info: list given by function get.M3.var().
##   llx, urx: Lower and upper x-coordinate bounds for the subsetted
##     grid in units appropriate to the model projection.  Defaults
##     are the current boundaries of the x range.
##   lly, ury: Lower and upper y-coordinate bounds for the subsetted
##     grid in units appropriate to the model projection.  Defaults
##     are the current boundaries of the y range.
##   ldatetime, udatetime: Starting and ending date-times (either Date
##     or POSIX class).  Defaults are the current boundaries of the
##     date-time range.
##   hz.strict: If TRUE (default), to be allowed in the subset,
##     the whole grid cell must fit within the bounds given by llx, urx,
##     lly, and ury.  If FALSE, grid cells will be included in the
##     subset if any portion of the grid cell's area falls within the
##     given bounds.
##
## RETURNS: Subsetted list, with appropriate elements altered to
##   reflect the subsetting that has taken place.
##
##
## NOTE: If the user wants to subset the variable by row, column,
##   layer, or time step number, this can be accomplished easily
##   using standard R methods for subsetting the array of variable
##   values.  This function was written to help the user who does not
##   know the row, column, or time step numbers, but who wants to subset
##   according to human-readable dates and times or according to
##   projection units.
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-05-20
## ###########################################################
var.subset <- function(var.info, llx, urx, lly, ury,
                       ldatetime, udatetime, hz.strict=TRUE){

  ## #############
  ## DEAL WITH COLUMNS FIRST.

  ## How many columns are there?
  num.columns <- length(var.info$x.cell.ctr)

  if (num.columns > 1){

    ## Want to find the left and right bounds (cell edges) for each
    ## column.  For this, we need the cell width.
    cell.width <- var.info$x.cell.ctr[2] - var.info$x.cell.ctr[1]
    lbd <- var.info$x.cell.ctr - (cell.width/2.0)
    ubd <- var.info$x.cell.ctr + (cell.width/2.0)
    rm(cell.width)

    ## If llx/urx is missing, set to preserve current boundaries.
    if (missing(llx))
      llx <- min(lbd)
    if (missing(urx))
      urx <- max(ubd)

    ## Check to make sure left limit is not greater than right limit.
    if (llx > urx)
      stop("Lower limit in the x direction is greater than upper limit.")

    ## If hz.strict=TRUE, find the columns for which both the
    ## left and right sides fit inside the specified x range. If
    ## hz.strict=FALSE, then only some portion of the column needs
    ## to fall within the bounds.
    if (hz.strict)
      which.columns <- which( (lbd >= llx) & (ubd <= urx) )
    else
      which.columns <- which( (ubd > llx) & (lbd < urx) )
  }

  else{
    which.columns <- 1
    if ( (!missing(llx)) || (!missing(urx)) )
      message("Only one column for this variable, llx and urx input ignored.")
  }
  ## #############


  ## #############
  ## DEAL WITH ROWS.

  ## How many rows are there?
  num.rows <- length(var.info$y.cell.ctr)

  if (num.rows > 1){

    ## Want to find the bottom and top bounds (cell edges) for each
    ## row.  For this, we need the cell width.
    cell.width <- var.info$y.cell.ctr[2] - var.info$y.cell.ctr[1]
    lbd <- var.info$y.cell.ctr - (cell.width/2.0)
    ubd <- var.info$y.cell.ctr + (cell.width/2.0)
    rm(cell.width)

    ## If lly/ury is missing, set to preserve current boundaries.
    if (missing(lly))
      lly <- min(lbd)
    if (missing(ury))
      ury <- max(ubd)

    ## Check to make sure left limit is not greater than right limit.
    if (lly > ury)
      stop("Lower limit in the y direction is greater than upper limit.")
    
    ## If hz.strict=TRUE, find the rows for which both the
    ## bottom and top sides fit inside the specified y range. If
    ## hz.strict=FALSE, then only some portion of the row needs
    ## to fall within the bounds.
    if (hz.strict)
      which.rows <- which( (lbd >= lly) & (ubd <= ury) )
    else
      which.rows <- which( (ubd > lly) & (lbd < ury) )
  }

  else{
    which.rows <- 1
    if ( (!missing(lly)) || (!missing(ury)) )
      message("Only one row for this variable, lly and ury input ignored.")
  }
  ## #############


  ## #############
  ## DEAL WITH DATE-TIMES.

  ## How many date-time steps are there?
  num.datetimes <- length(var.info$datetime)
  ## If this length is 0, then var.info$datetime is NULL.  This means
  ## that the file is time-independent, and we cannot subset the
  ## date-time steps.
  if (num.datetimes == 0)
    time.indep <- TRUE
  else
    time.indep <- FALSE


  if (num.datetimes > 1){

    ## If date-time limits are missing, set to preserve current boundaries.
    if (missing(ldatetime))
      ldatetime <- min(var.info$datetime)
    if (missing(udatetime))
      udatetime <- max(var.info$datetime)


    ## Check to see if the date-time limits are in Date format.  If so,
    ## make them into a POSIX format date.  For the lower limit, this
    ## would mean a time stamp at midnight (beginning of the given
    ## date).  For the upper limit, this would mean a time stamp at 23:59:59
    ## (last part of the given date).
    if ("Date" %in% class(ldatetime))
      ldatetime <- combine.date.and.time(date=ldatetime, time="00:00:00")
    if ("Date" %in% class(udatetime))
      udatetime <- combine.date.and.time(date=udatetime, time="23:59:59")


    ## Check to make sure upper limit is not less than lower limit.
    if (ldatetime > udatetime)
      stop("Lower date-time limit is greater than upper date-time limit.")


    ## Find the columns for which both the left and right sides fit
    ## inside the specified x range.
    which.datetimes <- which( (ldatetime <= var.info$datetime)
                             & (var.info$datetime <= udatetime) )
  }

  
  else if (!time.indep){
    which.datetimes <- 1
    if ( (!missing(ldatetime)) || (!missing(udatetime)) )
      message("Only one date-time step for this variable, ldatetime and udatetime input ignored.")
  }


  else{
    if ( (!missing(ldatetime)) || (!missing(udatetime)) )
      message("Variable is time-independent, ldatetime and udatetime input ignored.")
  }
  ## #############


  ## #############
  ## SUBSET THE VARIOUS ELEMENTS OF THE INPUT VARIABLE INFO LIST, AND
  ## RETURN THE SUBSETTED INFO.
  
  ## Subset the array first, recognizing that time-independent files
  ## cannot be subsetted on the basis of date-time.
  if (!time.indep){
    subset.data <- var.info$data[which.columns, which.rows, , which.datetimes]
    dim(subset.data) <- c(length(which.columns), length(which.rows), dim(var.info$data)[3], length(which.datetimes))
  }
  else{
    subset.data <- var.info$data[which.columns, which.rows, ]
    dim(subset.data) <- c(length(which.columns), length(which.rows), dim(var.info$data)[3])
  }


  ## Subset the x.cell.ctr and y.cell.ctr list elements.
  subset.x.cell.ctr <- var.info$x.cell.ctr[which.columns]
  subset.y.cell.ctr <- var.info$y.cell.ctr[which.rows]
  ## Subset the rows and columns list elements.
  subset.cols <- var.info$cols[which.columns]
  subset.rows <- var.info$rows[which.rows]
  ## If not time-independent, subset the datetime element.
  if (!time.indep)
    subset.datetime <- var.info$datetime[which.datetimes]

  
  ## Form the list for the subsetted variable info.
  if (!time.indep)
    subset.list <- list(data=subset.data, data.units=var.info$data.units,
                        x.cell.ctr=subset.x.cell.ctr,
                        y.cell.ctr=subset.y.cell.ctr,
                        hz.units=var.info$hz.units,
                        rows=subset.rows, cols=subset.cols,
                        layers=var.info$layers,
                        datetime=subset.datetime)
  else
    subset.list <- list(data=subset.data, data.units=var.info$data.units,
                        x.cell.ctr=subset.x.cell.ctr,
                        y.cell.ctr=subset.y.cell.ctr,
                        hz.units=var.info$hz.units,
                        rows=subset.rows, cols=subset.cols,
                        layers=var.info$layers)
  ## #############


  ## Return subsetted variable info to user.
  return(subset.list)
}
