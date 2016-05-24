#' Time Projection
#'
#' Project dates to lower dimensional subspace.
#' Extracts components year, month, yday, mday, hour, minute, weekday,
#' bizday and season from a date object
#'
#' @param dates date or datetime objects
#' @param size either "narrow" or "wide".  If narrow, returns a data frame
#'    containing the projections as column variables using factors.
#'    If wide, returns a sparse matrix containing the projections as column
#'    variables using 0-1 variables
#' @param holidays argument to determine which days are considered holidays,
#'    affecting the business day projection.  By default uses holidayNYSE()
#'    provided by the timeDate package, or can be specified as a vector of
#'    strings representing dates in the yyyy-mm-dd format
#' @param as.numeric logical only used when size = "narrow".  Returns the
#'    columns as numeric values instead of factors
#' @param drop logical.  If true, drop any column that only has 1 level or only
#'    1 unique element in it
#' @examples
#'    dates = timeSequence(from = "2001-01-01", to = "2004-01-01", by = "day")
#'    projectDate(as.Date(dates))
#' @export
projectDate = function(dates, size = c("narrow", "wide"),
                       holidays = holidayNYSE(year = unique(year(dates))),
                       as.numeric = F, drop = T) {
    size = match.arg(size)
    year = year(dates)
    month = month(dates)
    yday = yday(dates)
    mday = mday(dates)
    yweek = floor((yday - 1) / 7) + 1
    mweek = floor((mday - 1) / 7) + 1
    hour = hour(dates)
    minute = minute(dates)

    if (!as.numeric | size == "wide") {
        year = factor(year, ordered = T)
        month = factor(month, levels = 1:12, ordered = T)
        yday = factor(yday, levels = 1:366, ordered = T)
        mday = factor(mday, levels = 1:31, ordered = T)
        yweek = factor(yweek, levels = 1:53, ordered = T)
        mweek = factor(mweek, levels = 1:5, ordered = T)
        hour = factor(hour, levels = 0:23, ordered = T)
        minute = factor(minute, levels = 0:59, ordered = T)
    }
    
    weekday = factor(weekdays(dates), levels = c("Sunday", "Monday", "Tuesday", "Wednesday",
                                                 "Thursday", "Friday", "Saturday"),
                     ordered = T)
    bizday = factor(is.Bizday(dates, holidays))
    season = factor(getSeason(dates), levels = c("Winter", "Spring", "Summer", "Fall"),
                    ordered = T)
    raw = data.frame(year = year,
                     month = month,
                     yweek = yweek,
                     mweek = mweek,
                     yday = yday,
                     mday = mday,
                     hour = hour,
                     minute = minute,
                     weekday = weekday,
                     bizday = bizday,
                     season = season)
    if (drop) {
        redundantCols = rep(F, ncol(raw))
        for (i in 1:ncol(raw)) {
            if (all(diff(as.numeric(raw[,i])) == 0)) redundantCols[i] = T
        }
        #redundantCols = apply(raw, 2, function(j) { nlevels(j) == 1 })
        raw = raw[,!redundantCols]
    }
    if (size == "narrow") return (raw)
    if (size == "wide") {
        return (sparse.model.matrix(~ . -1, raw))
    }
}

# Source: http://stackoverflow.com/questions/9500114/find-which-season-a-particular-date-belongs-to
getSeason <- function(dates) {
    WS <- as.Date("2012-12-15", format = "%Y-%m-%d") # Winter Solstice
    SE <- as.Date("2012-3-15",  format = "%Y-%m-%d") # Spring Equinox
    SS <- as.Date("2012-6-15",  format = "%Y-%m-%d") # Summer Solstice
    FE <- as.Date("2012-9-15",  format = "%Y-%m-%d") # Fall Equinox

    # Convert dates from any year to 2012 dates
    d <- as.Date(strftime(dates, format="2012-%m-%d"))

    ifelse (d >= WS | d < SE, "Winter",
      ifelse (d >= SE & d < SS, "Spring",
        ifelse (d >= SS & d < FE, "Summer", "Fall")))
}

is.Bizday = function(x, holidays) 
{
    char.x = substr(as.character(x), 1, 10)
    char.h = substr(as.character(holidays), 1, 10)
    Weekday = as.integer(isWeekday(x, wday = 1:5))
    nonHoliday = as.integer(!(char.x %in% char.h))
    bizdays = as.logical(Weekday * nonHoliday)
    return (bizdays)
}

#' @method weekdays timeDate
#' @S3method weekdays timeDate
weekdays.timeDate = function(x, abbreviate=FALSE) {
    weekdays(as.Date(x), abbreviate)
}
