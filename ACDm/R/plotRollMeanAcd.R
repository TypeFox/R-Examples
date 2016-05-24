plotRollMeanAcd <- function(durations, window = 500){
  
  rollingMeans <- group <- NULL
  
  dur <- durations$durations
  time <- durations$time[(window):length(dur)]
  if(!("POSIXlt" %in% class(time))) time <- as.POSIXlt(time)
    
  df <- data.frame(time = time, rollingMeans = zoo::rollmean(dur,window), group = time$year + time$yday/365)
  
  g <- ggplot(df, aes(x = time,y = rollingMeans, group = group))
  g + geom_line()
}