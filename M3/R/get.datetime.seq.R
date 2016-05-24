## ###########################################################
## PURPOSE: Read the date-time steps in the Models3-formatted file.
##   Put these into datetime format.
##
## INPUT:
##   file: File name of Models3-formatted file of interest.
##
## RETURNS: List of datetimes included in the Models3-formatted file.
##
## ASSUMES: This code assumes that the time step is not negative.
##   (For instance, the Models3 I/OAPI does allow for negative time
##   steps, but these negative time steps will NOT be handled properly
##   by this function.)
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-06-02
## ###########################################################
get.datetime.seq <- function(file){

  ## Open netCDF file which has the projection we want to use.
  nc <- nc_open(file)


  ## Get information about the time step increment.
  tstep.incr <- ncatt_get(nc, varid=0, attname="TSTEP")$value

  ## Test to see whether the file is time-independent (in which case,
  ## the time increment will be 0).
  if (tstep.incr==0){
    ## Close the Models3 file, issue warning, and exit the function.
    nc <- nc_close(nc)
    warning("Time step increment is zero.  This appears to be a time-independent file.")
    return(NULL)
  }

  
  ## If the time increment is not zero, then we assume there are time
  ## steps for the variables in the file.

  ## Find the starting date and time.
  M3.start.date <- ncatt_get(nc, varid=0, attname="SDATE")$value
  start.date <- decipher.M3.date(M3.start.date)
  M3.start.time <- ncatt_get(nc, varid=0, attname="STIME")$value
  start.time <- decipher.M3.time(M3.start.time)

  ## Combine the starting date and time to form datetime object (POSIX
  ## class).
  start.datetime <- combine.date.and.time(date=start.date, time=start.time)

  
  ## Find the increment separating the time steps.
  tstep.incr.list <- decipher.M3.time(tstep.incr)
  ## To find the increment of the time step in seconds.
  tstep.in.secs <- (tstep.incr.list$hrs*60*60) + (tstep.incr.list$mins*60) + tstep.incr.list$secs

  
  ## How many datetimes are there?  (What is the length of the time
  ## dimension?)
  num.time.steps <- nc$dim$TSTEP$len

  
  ## Now get a sequence.
  datetime.seq <- seq.POSIXt(from=start.datetime, by=tstep.in.secs,
                             length.out=num.time.steps)

  
  ## Close the Models3 file.
  nc <- nc_close(nc)
  rm(nc)


  ## Return the date-time sequence we've developed.
  return(datetime.seq)
}
