## ###########################################################
## PURPOSE: Combine date and time to obtain date-time in POSIX format.
##
## INPUTS:
##   date: Date in Date format or as character string in format "YYYY-MM-DD".
##   time: Time as list with hrs, mins, and secs components or as
##      character string in "HH:MM:SS" (with hours 00-23).
##
## RETURNS: A date-time in POSIX format.
##
## ASSUMES: This code assumes that the time is not negative.  (For
##   instance, the Models3 I/OAPI does allow for negative time steps,
##   but these negative time steps will NOT be handled properly by
##   this function.)
##
##
## REVISION HISTORY:
##   Original release: Jenise Swall, 2011-05-19
## ###########################################################
combine.date.and.time <- function(date, time){

  ## Check whether time is a list like that returned by the
  ## decipher.M3.time() function.
  if (is.list(time))  
    datetime <- strptime(paste(as.character(date), " ", time$hrs, ":",
                               time$mins, ":", time$secs, sep=""),
                         format="%Y-%m-%d %H:%M:%S", tz="GMT")

  ## Otherwise, assume time is a character string of form HH:MM:SS.
  else
    datetime <- strptime(paste(as.character(date), " ", time, sep=""),
                         format="%Y-%m-%d %H:%M:%S", tz="GMT")

  return(datetime)
}
