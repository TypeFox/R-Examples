## ###########################################################
## PURPOSE: Decipher Models3 time format (HHMMSS) into hours, minutes,
##   and seconds.
##
## INPUT: Models3 time (numeric) in the format HHMMSS
##
## RETURNS: List with hrs (hours), mins (minutes), and secs (seconds)
##   components.
##
## ASSUMES: This code assumes that the time is not negative.  (For
##   instance, the Models3 I/OAPI does allow for negative time steps,
##   but these negative time steps will NOT be handled properly by
##   this function.)
##
## NOTE: The Models3 time is an integer, so we can't just extract the
##   first 2 characters, next 2 characters, etc.  If the time step is
##   one hour, then the time we extract will be 100, not 000100.
##
##
## RELEASE HISTORY:
##   Original release: Jenise Swall, 2011-05-19
## ###########################################################
decipher.M3.time <- function(M3.time){

  ## Find number of hours.
  hrs <- trunc(M3.time/10000)

  ## Find number of minutes.
  mins <- trunc( (M3.time - (hrs * 10000)) / 100 )

  ## Find number of seconds
  secs <- M3.time - (hrs*10000) - (mins * 100)

  return(list(hrs=hrs, mins=mins, secs=secs))
}
