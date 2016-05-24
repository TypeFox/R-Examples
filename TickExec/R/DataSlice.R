#### for given data.frame with a column named 'time' ####
#### returns a sub data.frame between time1 and min(time2, time1 + last) ####

DataSlice <- function(df, time1 = 000001, time2 = 235959, last = 24 * 3600) {
  ## formalize argument ##
  if (class(df) != 'data.frame') {
    stop('Need dataframe as input.')
  } 
  
  colnames(df) <- toupper(colnames(df))
  if (is.na(match('TIME', colnames(df))) == TRUE) {
    stop ('No time column in dataframe.')
  }
  
  ## construct time interval ##
  deltaTime = min(TimeDiff(time1, time2), last)
  endTime = TimeAdd(time1, deltaTime)
  
  idx <- as.numeric(df$TIME)
  out <- df[which(idx >= time1 & idx <= endTime), ]
  
  if (nrow(out) < 1) {
    return (NA)
  } else {
    return (out)
  }
}