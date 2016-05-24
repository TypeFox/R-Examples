plotHistAcd <- function(durations, windowunit = "mins", window = 1){
  
  time <- NULL
  
  if(!("POSIXlt" %in% class(durations$time))) durations$time <- as.POSIXlt(durations$time)
  windowunit <- match.arg(windowunit, c("secs", "mins", "hours", "days"))
  
  timefactor <- switch(windowunit,
                       secs = trunc(durations$time, units = windowunit) - durations$time$sec %% window,
                       mins = trunc(durations$time, units = windowunit) - (durations$time$min %% window) * 60,
                       hours = trunc(durations$time, units = windowunit) - (durations$time$hour %% window) * 3600,
                       days = trunc(durations$time, units = windowunit) - (durations$time$yday %% window) * 86400)
  
  meandur <- tapply(durations$durations, timefactor, mean)
  
  df <- data.frame(time = as.POSIXlt(names(meandur)), meandur = meandur)
  
  g <- ggplot(df, aes(x = time, y = meandur))
  g + geom_bar(stat = "identity") + ggtitle(paste(window, windowunit, "means"))
}