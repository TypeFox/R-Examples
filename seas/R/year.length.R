"year.length" <-
function(x, calendar) {
  if (inherits(x, c("POSIXct", "Date")))
    x <- format(x, "%Y")
  x <- as.integer(x)
  if (missing(calendar))
    calendar <- NULL
  if (!is.null(calendar)) {
    if (!calendar %in% c("360_day", "365_day", "366_day", "julian",
                         "noleap", "all_leap"))
      stop("Invalid 'calendar' argument;\nmust be 'NA' (for Gregorian), 'julian', '360_day', '365_day' (no leap), or '366_day' (all leap)")
    if (calendar=="360_day")
      return(rep(360, length(x)))
    if (calendar %in% c("365_day", "noleap"))
      return(rep(365, length(x)))
    if (calendar %in% c("366_day", "all_leap"))
      return(rep(366, length(x)))
    if (calendar == "julian")
      return(rep(365, length(x)) + ifelse(x %% 4 == 0, 1, 0))
  }  # default is proleptic Gregorian
  return(rep(365, length(x)) + ifelse(x %% 4 == 0 & x %% 100 != 0 | x %% 400 == 0, 1, 0))
}
