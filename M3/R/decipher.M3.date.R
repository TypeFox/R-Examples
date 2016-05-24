## ###########################################################
## PURPOSE: Decipher Models3 date format (YYYYDDD) into R's Date
##   class.
##
## INPUT: Models3 date (numeric) in the format YYYYDDD, where DDD is
##   a Julian day (since the beginning of YYYY).
##
## RETURNS: Starting date YYYYDDD in R's Date class.
##
##
## NOTE: The Models3 date is an integer, so we can't just extract the
##   first 4 characters, next 3 characters, etc.
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-05-19
## ###########################################################
decipher.M3.date <- function(M3.date){

  ## Find the year.
  yr <- trunc(M3.date/1000)

  ## Find number of minutes.
  julian.day <- M3.date %% 1000

  ## The first day of this year is our base date.
  jan1.yr <- as.Date(paste(yr, "-01-01", sep=""), format="%Y-%m-%d")

  ## Pass to as.Date function the number of days since the "origin"
  ## date.  Our origin date is Jan. 1, YYYY, which is stored as a Date
  ## object in jan1.yr.  Note that Jan. 1, YYYY is Julian day 001.
  my.date <- as.Date(julian.day-1, origin=jan1.yr) 

  return(my.date)
}
