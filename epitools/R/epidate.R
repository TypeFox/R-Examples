##in projgress
"epidate" <-
  function(x, format = "%m/%d/%Y", cal.dates = FALSE,
           before = 7, after = 7, sunday = TRUE){
    dates  <- as.Date(x, format = format)
    julian <- julian(dates)
    posixlt <- as.POSIXlt(dates)
    mday    <- posixlt$mday ##1-31: day of the month
    mon     <- posixlt$mon ##0-11: months after the first of the year.
    month   <- months(dates) ##January, February, ....
    month2  <- substr(month, 1, 3) ##Jan, Feb, ...
    if(sunday){
      week <- format(dates, format = "%U") ##0 1 2 3 ... 51 52 53
      firstday <- "Sunday"
    } else {
      week <- format(dates, format = "%W")
      firstday <- "Monday"
    }
    names(firstday) <- "1st day of week for $week number:"
    year    <- posixlt$year + 1900 ##Years0
    yr      <- substr(as.character(year), 3, 4)
    wday    <- posixlt$wday ##0-6 day of the week, starting on Sunday.
    weekday <- weekdays(dates) ##Sunday Monday Tuesday
    wkday   <- substr(weekday, 1, 3) ##Sun Mon Tue
    yday    <- posixlt$yday ##0-365: day of the year
    quarter <- quarters(dates)
    if(cal.dates==TRUE){
      cdates  <- seq(from = min(dates, na.rm = TRUE) - before,
                     to = max(dates, na.rm = TRUE) + after, by = 1)
      cjulian  <- julian(cdates)
    } else {
      cdates <- "Not reported: To report, set cal.dates=TRUE"
      cjulian <- "Not reported: To report, set cal.dates=TRUE"
    }
    list(dates  = dates,
         julian = julian,
         mday    = mday,
         mon     = mon,
         month   = month,
         month2  = month2,
         firstday = firstday,       
         week    = week,
         year    = year,
         yr      = yr,
         wday    = wday,
         weekday = weekday,
         wkday   = wkday,
         yday    = yday,
         quarter = quarter,
         cdates  = cdates,
         cjulian  = cjulian
         )
  }
