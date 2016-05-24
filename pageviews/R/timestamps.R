#'@title Convert time objects to function with pageviews functions
#'@description \code{pageview_timestamps} converts \code{\link{Date}} and \code{\link{POSIXlt}} and ct
#'objects to work nicely with the \code{start} and \code{end} parameters in pageviews functions.
#'
#'@param timestamps a vector of Date, POSIXlt or POSIXct objects.
#'
#'@param first whether to, if \code{timestamps} is of date objects, assume the
#'first hour in a day (TRUE) or the last (FALSE). TRUE by default.
#'
#'@return a character vector containing timestamps that can be used with \code{\link{article_pageviews}} et al.
#'
#'@seealso \code{\link{article_pageviews}} and \code{\link{project_pageviews}}, where you
#'can make use of this function.
#'
#'@examples
#'# Using a Date
#'pageview_timestamps(Sys.Date())
#'
#'# Using a POSIXct object
#'pageview_timestamps(Sys.time())
#'
#'@export
pageview_timestamps <- function(timestamps = Sys.Date(), first = TRUE) {
  
  # Check type
  classes <- c("Date", "POSIXlt", "POSIXct")
  if(!any(class(timestamps) %in% classes)){
    stop("'timestamps' must be of type Date, POSIXlt or POSIXct")
  }
  is_date <- ifelse(("Date" %in% class(timestamps)), TRUE, FALSE)
  
  # Convert to character and strip out unnecessary chars
  timestamps <- gsub(x = timestamps, pattern = "( |-|:)", replacement = "")
  
  # If it's a date, it'll only be 6 chars long, so we have to add 4 more
  # using "first" to determine what they should be.
  if(is_date){
    if(first){
      return(paste0(timestamps, "0100"))
    }
    return(paste0(timestamps, "2400"))
  }
  
  # Otherwise it was a POSIX object so we have the opposite problem.
  timestamps <- substring(timestamps, 0, 10)
  return(paste0(timestamps, "00"))
}
