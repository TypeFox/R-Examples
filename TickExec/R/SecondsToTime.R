#### return numeric time of the day ####

SecondsToTime <- function(seconds) {
  ## formalize argument ##
  if (seconds > 24*3600) {
    stop('Can not cope more than one day.')
  }
  
  h = seconds %/% 3600
  ms = seconds %% 3600
  m = ms %/% 60
  s = ms %% 60
  
  out <- sprintf('%02d%02d%02d', h, m, s)
  
  return(as.numeric(out))
}