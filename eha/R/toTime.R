toTime <- function(dates){
    if (is.numeric(dates)) return(dates)
    dates <- as.character(dates)
    c(as.Date(dates) - as.Date("0000-01-01") ) / 365.2425
  }
