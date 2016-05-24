time365toDate <- function(x)
{
  year  <- trunc(x)
  md <- .time365md[round(1 + 365*(x%%1)), ]
  return(as.Date(paste(year, md$month, md$day, sep = '-')))
}
