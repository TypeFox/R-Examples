#' Convert MMWRweek to Date
#' 
#' Computes the Date from the MMWR year, week, and day.
#' 
#' @param MMWRyear numeric vector of years
#' @param MMWRweek numeric vector of weeks
#' @param MMWRday numeric vector of days, defaults to a vector of 1s
#' @return Date vector of dates associated with MMWR year, week, and day
#' @author Jarad Niemi \email{niemi@@iastate.edu}
#' @seealso \code{\link{MMWRweek}}
#' @export
#' @examples
#' MMWRweek2Date(MMWRyear=2015,MMWRweek=36,MMWRday=3)
MMWRweek2Date = function(MMWRyear, MMWRweek, MMWRday=NULL) {
  stopifnot(all(is.numeric(MMWRyear)))
  stopifnot(all(is.numeric(MMWRweek)))
  stopifnot(all(0 < MMWRweek & MMWRweek < 54))

  stopifnot(length(MMWRyear) == length(MMWRweek))
  
  if (is.null(MMWRday)) MMWRday = rep(1,length(MMWRweek))
  stopifnot(all(0 < MMWRday & MMWRday < 8))
  
  jan1 = start_date(MMWRyear)
  return(jan1 + (MMWRweek-1)*7 + MMWRday-1)
}
