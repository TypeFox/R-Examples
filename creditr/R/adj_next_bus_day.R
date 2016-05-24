#' Adjust to next business day.
#' 
#' \code{adj_next_bus_day} gets the next business day following 5D bus day
#' convention.
#' 
#' @param date, a \code{Date} type.
#'   
#' @return Date adjusted to the following business day

adj_next_bus_day <- function(date){
  
  for(i in 1:length(date)){
    
    dateWday <- as.POSIXlt(as.Date(date[i]))$wday
    
    ## change date to the most recent weekday if necessary
    
    if (dateWday == 0){
      date[i] <- date[i] + 1
    } else if (dateWday == 6) {
      date[i] <- date[i] + 2
    }
  }
  return(date)
}
