
date2day <- function(dates, start=NULL, tz="", ...)
{
  dates <- as.POSIXlt(dates, tz=tz, ...)
  if (is.null(start))
    start <- min(dates)
  else
    start <- as.POSIXlt(start, tz=tz, ...)
  out <- as.numeric(difftime(dates, start, units="days"))
  return(out)
}
