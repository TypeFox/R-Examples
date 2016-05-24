# this function was copied from the package date April 2014
# we copied the function to our package, to reduce the number 
# of packages on which birdring is dependent.
#--------------------------------------------------------
# For other functions to handle date variables, see
# Terry Therneau and Thomas Lumley and Kjetil Halvorsen
# and Kurt Hornik (2012). date: Functions for handling
# dates. R package version 1.2-33.
#-----------------------------------------------------------

mdy.date <- function (month, day, year, nineteen = TRUE, fillday = FALSE, 
          fillmonth = FALSE) 
{
  temp <- any((month != trunc(month)) | (day != trunc(day)) | 
                (year != trunc(year)))
  if (!is.na(temp) && temp) {
    warning("Non integer input values were truncated in mdy.date")
    month <- trunc(month)
    day <- trunc(day)
    year <- trunc(year)
  }
  if (nineteen) 
    year <- ifelse(year < 100, year + 1900, year)
  temp <- numeric(length(month + day + year))
  month <- month + temp
  day <- day + temp
  year <- year + temp
  if (fillmonth) {
    temp <- is.na(month)
    month[temp] <- 7
    day[temp] <- 1
  }
  if (fillday) 
    day[is.na(day)] <- 15
  month[month < 1 | month > 12] <- NA
  day[day < 1] <- NA
  year[year == 0] <- NA
  year <- ifelse(year < 0, year + 1, year)
  tyear <- ifelse(month > 2, year, year - 1)
  tmon <- ifelse(month > 2, month + 1, month + 13)
  julian <- trunc(365.25 * tyear) + trunc(30.6001 * tmon) + 
    day - 715940
  temp <- trunc(0.01 * tyear)
  save <- ifelse(julian >= -137774, julian + 2 + trunc(0.25 * 
                                                         temp) - temp, julian)
  year <- ifelse(month == 12, year + 1, year)
  month <- ifelse(month == 12, 1, month + 1)
  day <- 1
  tyear <- ifelse(month > 2, year, year - 1)
  tmon <- ifelse(month > 2, month + 1, month + 13)
  julian <- trunc(365.25 * tyear) + trunc(30.6001 * tmon) + 
    day - 715940
  temp <- trunc(0.01 * tyear)
  save2 <- ifelse(julian >= -137774, julian + 2 + trunc(0.25 * 
                                                          temp) - temp, julian)
  temp <- as.integer(ifelse(save2 > save, save, NA))
  attr(temp, "class") <- "date"
  temp
}