#' Find start date for a calendar year
#' 
#' Finds the state date given a numeric calendar year
#' 
#' @param year integer vector of four digit years
#' @return Date vector for January 1st of the calendar year
#' @seealso \code{\link{get_start_date}}
#' @author Jarad Niemi \email{niemi@@iastate.edu}
start_date = function(year) {
  # Finds start state for this calendar year
  jan1 = as.Date(paste(year, '-01-01', sep=''))
  wday = as.numeric(MMWRweekday(jan1))
  jan1 - (wday-1) + 7*(wday>4)
}

#' Finds the start date for the year associated with date.
#' 
#' @param date vector which can be coerced to class \code{Date}
#' @return Date vector for start date of MMWR year associated with date
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @seealso \code{\link{start_date}}
get_start_date <- function(date) {
  year = as.integer(format(as.Date(date),'%Y'))
  
  # Find start date for current and previous calendar years
  sd_prev    = start_date(year-1)
  sd_current = start_date(year  )
  sd_next    = start_date(year+1)
  
  before = date <  sd_current
  after  = date >= sd_next
  
  out = sd_current
  out[before] = sd_prev[before]
  out[after]  = sd_next[after]
  return(out)
}

