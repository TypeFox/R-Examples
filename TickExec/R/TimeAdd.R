#### return time2 in numeric ####

TimeAdd <- function(time1, increase) {
  ## formalize argument ##
  time1 <- as.numeric(time1)
  
  ## calculate end time ##
  seconds1 = TimeDiff(0, time1)
  seconds2 = seconds1 + increase
  if (seconds2 > 24 * 3600) {
    stop('Can not cope more than 1 day.')
  }
  time2 <- SecondsToTime(seconds2)
  
  return(time2)
}