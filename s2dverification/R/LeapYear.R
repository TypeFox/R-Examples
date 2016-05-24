LeapYear <- function(year) {
  leap <- FALSE
  if (year %% 4 == 0) {
    leap <- TRUE
    if (year %% 100 == 0) {
      leap <- FALSE
      if (year %% 400 == 0) {
        leap <- TRUE
      }
    } 
  }
  
  #
  # Output
  #
  leap
}
