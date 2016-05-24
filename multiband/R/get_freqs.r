get_time_range <- function(tms) {
  time_range <- 0
  for(tm in tms){
    time_range <- max(max(tm[,1]) - min(tm[,1]),time_range)
  }
  return(time_range)
}

get_freq_del <- function(tms) {
  time_range <- get_time_range(tms)
  return(.1/time_range)
}

#' Construct grid of frequencies
#' 
#' \code{get_freqs} constructs a grid of frequencies from a specified minimum and maximum period and frequency grid length.
#' 
#' @param period_min minimum period
#' @param period_max maximum period
#' @param freq_del gride size in frequency
#' @export
get_freqs <- function(period_min=1,period_max=100,freq_del=.1/4000) {
  freq_max <- 1 / period_min
  freq_min <- 1 / period_max
  ## Richards "On Machine . . ." recommends grid
  ## with intervals of length .1/(max(time) - min(time))
  ## (max(time) - min(time)) ~ 4000 for most OGLE curves
  ## construct sequence and convert frequency to radians
  return(2*pi*seq(freq_min,freq_max,freq_del))
}
