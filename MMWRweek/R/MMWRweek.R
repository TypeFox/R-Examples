#' MMWR day, week, and year 
#'
#' This function returns the MMWR day, week, and year for the Date(s) provided.
#' 
#' The first day of any MMWR week is Sunday. MMWR week numbering is sequential 
#' beginning with 1 and incrementing with each week to a maximum of 52 or 53. 
#' MMWR week #1 of an MMWR year is the first week of the year that has at least 
#' four days in the calendar year. For example, if January 1 occurs on a Sunday, 
#' Monday, Tuesday or Wednesday, the calendar week that includes January 1 would 
#' be MMWR week #1. If January 1 occurs on a Thursday, Friday, or Saturday, the 
#' calendar week that includes January 1 would be the last MMWR week of the previous 
#' year (#52 or #53). Because of this rule, December 29, 30, and 31 could potentially 
#' fall into MMWR week #1 of the following MMWR year.
#' 
#' @param date vector which can be coerced to class \code{Date}
#' @return data.frame with elements MMWRday (of the week), MMWRweek, and MMWRyear
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @seealso \code{\link{MMWRweekday}}, \code{\link{MMWRweek2Date}}
#' @references \url{http://wwwn.cdc.gov/nndss/document/MMWR_Week_overview.pdf}
#' @export 
#' @examples
#' y = as.Date(paste(1999:2011, "-12-31", sep = ""))
#' cbind(y, MMWRweek(y))
MMWRweek = function(date) {
  date = as.Date(date)
  #year = as.numeric(format(date, '%Y'))
  start_date = get_start_date(date)

  return(data.frame(MMWRyear = as.numeric(format(start_date+7, '%Y')),
                    MMWRweek = round(as.numeric(date - start_date - 3) / 7)+1, # why 3? 
                    MMWRday  = as.numeric(MMWRweekday(date)))
  )
} 
