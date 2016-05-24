#' Scale event counts
#'
#' Scale event counts based on the unit of analysis.
#' 
#' @aliases NormEventCounts
#' @param x data.frame, a GDELT data.frame.
#' @param unit.analysis character, default is country.day; other options: country.day, country.month, country.year, day, month, year
#' @param var.name character, base name for the new count variables
#' @return data.frame
#' @export
#' @details
#' For \code{unit.analysis}, day and country-day put out a data set where date
#' is of class \sQuote{date}.  All other options put out a data set where year
#' or month is integer (this needs to be unified in a later version).
#' @references
#' GDELT: Global Data on Events, Location and Tone, 1979-2012.  
#' Presented at the 2013 meeting of the International Studies Association
#' in San Francisco, CA.
#' \url{http://www.gdeltproject.org/}
#' @author 
#' \tabular{ll}{
#'   Oskar N.T. Thoms \tab \email{othoms@@princeton.edu}\cr
#'   Stephen R. Haptonstahl \tab \email{srh@@haptonstahl.org}\cr
#'   John Beieler \tab \email{jub270@@psu.edu}\cr
#' }
#' @examples
#' \dontrun{
#' GDELT.subset.data <- GetGDELT("2013-06-01", "2013-06-07", allow.wildcards=TRUE,
#'   filter=list(ActionGeo_CountryCode=c("AF", "US"), EventCode="14*"),
#'   local.folder="~/gdeltdata")
#' GDELT.normed.data <- NormEventCounts(x = GDELT.subset.data, 
#'   unit.analysis="day", 
#'   var.name="protest")}
NormEventCounts <- function(x, 
                            unit.analysis, 
                            var.name="norming_vars"){
  # files with total event counts are included in the package, but these could be stored on a server and downloaded as needed 
  if(missing(x)) stop("Must specify a data.frame")
  if(missing(unit.analysis)) stop("A unit of analysis must be specified")
  if(!(unit.analysis %in% c("country.day", 
                            "country.month", 
                            "country.year", 
                            "day", 
                            "month", 
                            "year"))) stop("Not a valid unit of analysis. Choose from: country.day, country.month, country.year, day, month, year")

  local.folder <- tempdir()

  if(unit.analysis == "country.day"){ # see here for code annotations
    ##The pattern for downloading repeats for each unit of analysis
    #Set destination folder
    dest <- paste(local.folder, "/norm_counts_daily_country.csv", sep="")
    #Download data from the GDELT server at UMN
    dl.daily.country.data <- download.file(url="http://data.gdeltproject.org/normfiles/daily_country.csv", 
                                        destfile=dest, 
                                        quiet=FALSE)
    #Read in the downloaded CSV
    daily.country.data <- read.csv(dest, col.names=c("day", "country", "total"))

    # code event counts
    x$count <- tapply(x$EventCode, 
                      paste(x$ActionGeo_CountryCode, x$SQLDATE), length)[paste(x$ActionGeo_CountryCode, x$SQLDATE)]
    # merge the two together
    x <- merge(x, daily.country.data, 
               by.x = c("ActionGeo_CountryCode", "SQLDATE"), 
               by.y = c("country", "day"), 
               all.x = TRUE)
    
    # code the normalized count as (count of selected events in the country-day) / (all events in the country-day)
    x$norm.count <- x$count/x$total
    
    # get the first and last date in the data
    range <- as.Date(as.character(range(x$SQLDATE)), format = "%Y%m%d")
    # create a vector with all dates in the range
    days <- format(seq(as.Date(range[1]), to = as.Date(range[2]), by = "1 day"), "%Y%m%d")
    
    # only keep the columns needed & subset to unique country-days
    x <- unique(subset(x, select = c("ActionGeo_CountryCode", "SQLDATE", "count", "norm.count")))
    # create complete time-series per country
    complete <- expand.grid(country = unique(daily.country.data$country), date = days)
    # create new names based on user input
    new.vars <- c(paste(var.name, ".count", sep = ""), paste(var.name, ".norm", sep = ""))
    # assign column names
    names(x) <- c("country", "date", new.vars)
    # merge data into complete time-series
    x <- merge(complete, x, by = c("country", "date"), all.x = TRUE)
    # code NAs as 0s
    x[is.na(x[, paste(var.name, ".norm", sep = "")]), new.vars] <- c(0, 0)
    # coerce into date format
    x[, "date"] <- as.Date(x[, "date"], "%Y%m%d") 
  } 
  if(unit.analysis == "country.month"){
    dest <- paste(local.folder, "/norm_counts_monthly_country.csv", sep="")
    dl.monthly.country.data <- download.file(url="http://data.gdeltproject.org/normfiles/monthly_country.csv", 
                                        destfile=dest, 
                                        quiet=FALSE)
    monthly.country.data <- read.csv(dest, col.names=c("month", "country", "total"))

    x$count <- tapply(x$EventCode, paste(x$ActionGeo_CountryCode, x$MonthYear), length)[paste(x$ActionGeo_CountryCode, x$MonthYear)]
    x <- merge(x, monthly.country.data, by.x = c("ActionGeo_CountryCode", "MonthYear"), by.y = c("country", "month"), all.x = TRUE)
    x$norm.count <- x$count/x$total
    range <- as.Date(as.character(range(x$SQLDATE)), format = "%Y%m%d") # get the first and last date in the data
    months <- as.integer(format(seq(as.Date(range[1]), to = as.Date(range[2]), by = "1 month"), "%Y%m"))
    x <- unique(subset(x, select = c("ActionGeo_CountryCode", "MonthYear", "count", "norm.count")))
    complete <- expand.grid(country = unique(monthly.country.data$country), month = months) 
    new.vars <- c(paste(var.name, ".count", sep = ""), paste(var.name, ".norm", sep = ""))
    names(x) <- c("country", "month", new.vars) 
    x <- merge(complete, x, by = c("country", "month"), all.x = TRUE)
    x[is.na(x[, paste(var.name, ".norm", sep = "")]), new.vars] <- c(0, 0)
  } 
  if(unit.analysis == "country.year"){ 
    dest <- paste(local.folder, "/norm_counts_yearly_country.csv", sep="")
    dl.yearly.country.data <- download.file(url="http://data.gdeltproject.org/normfiles/yearly_country.csv", 
                                        destfile=dest, 
                                        quiet=FALSE)
    yearly.country.data <- read.csv(dest, col.names=c("year", "country", "total"))

    x$count <- tapply(x$EventCode, paste(x$ActionGeo_CountryCode, x$Year), length)[paste(x$ActionGeo_CountryCode, x$Year)]
    x <- merge(x, yearly.country.data, by.x = c("ActionGeo_CountryCode", "Year"), by.y = c("country", "year"), all.x = TRUE)
    x$norm.count <- x$count/x$total
    range <- range(x$Year)
    x <- unique(subset(x, select = c("ActionGeo_CountryCode", "Year", "count", "norm.count")))
    complete <- expand.grid(country = unique(yearly.country.data$country), year = range[1]:range[2])
    new.vars <- c(paste(var.name, ".count", sep = ""), paste(var.name, ".norm", sep = ""))
    names(x) <- c("country", "year", new.vars)
    x <- merge(complete, x, by = c("country", "year"), all.x = TRUE) 
    x[is.na(x[, paste(var.name, ".norm", sep = "")]), new.vars] <- c(0, 0)
  } 
  if(unit.analysis == "day"){
    dest <- paste(local.folder, "/norm_counts_daily.csv", sep="")
    dl.daily.data <- download.file(url="http://data.gdeltproject.org/normfiles/daily.csv", 
                                        destfile=dest, 
                                        quiet=FALSE)
    daily.data <- read.csv(dest, col.names=c("day", "total"))

    x$count <- tapply(x$EventCode, x$SQLDATE, length)[as.factor(x$SQLDATE)]
    x <- merge(x, daily.data, by.x = "SQLDATE", by.y = "day", all.x = TRUE)		
    x$norm.count <- x$count/x$total
    range <- as.Date(as.character(range(x$SQLDATE)), format = "%Y%m%d")
    days <- format(seq(as.Date(range[1]), to = as.Date(range[2]), by = "1 day"), "%Y%m%d")
    x <- unique(subset(x, select = c("SQLDATE", "count", "norm.count")))
    complete <- expand.grid(date = days) 
    new.vars <- c(paste(var.name, ".count", sep = ""), paste(var.name, ".norm", sep = ""))
    names(x) <- c("date", new.vars)
    x <- merge(complete, x, by = "date", all.x = TRUE)
    x[is.na(x[, paste(var.name, ".norm", sep = "")]), new.vars] <- c(0, 0)
    x[, "date"] <- as.Date(x[, "date"], "%Y%m%d")
  } 
  if(unit.analysis == "month"){
    dest <- paste(local.folder, "/norm_counts_monthly.csv", sep="")
    dl.monthly.data <- download.file(url="http://data.gdeltproject.org/normfiles/monthly.csv", 
                                        destfile=dest, 
                                        quiet=FALSE)
    monthly.data <- read.csv(dest, col.names=c("month", "total"))

    x$count <- tapply(x$EventCode, x$SQLDATE, length)[as.factor(x$SQLDATE)]
    x <- merge(x, monthly.data, by.x = "MonthYear", by.y = "month", all.x = TRUE)
    x$norm.count <- x$count/x$total
    range <- as.Date(as.character(range(x$SQLDATE)), format = "%Y%m%d")
    months <- as.integer(format(seq(as.Date(range[1]), to = as.Date(range[2]), by = "1 month"), "%Y%m"))
    x <- unique(subset(x, select = c("MonthYear", "count", "norm.count")))
    complete <- expand.grid(month = months) 
    new.vars <- c(paste(var.name, ".count", sep = ""), paste(var.name, ".norm", sep = ""))
    names(x) <- c("month", new.vars)
    x <- merge(complete, x, by = "month", all.x = TRUE)
    x[is.na(x[, paste(var.name, ".norm", sep = "")]), new.vars] <- c(0, 0)
  } 
  if(unit.analysis == "year"){
    dest <- paste(local.folder, "/norm_counts_yearly.csv", sep="")
    dl.yearly.data <- download.file(url="http://data.gdeltproject.org/normfiles/yearly.csv", 
                                        destfile=dest, 
                                        quiet=FALSE)
    yearly.data <- read.csv(dest, col.names=c("year", "total"))

    x$count <- tapply(x$EventCode, x$Year, length)[as.factor(x$Year)]
    x <- merge(x, yearly.data, by.x = "Year", by.y = "year", all.x = TRUE)		
    x$norm.count <- x$count/x$total
    range <- range(x$Year)
    x <- unique(subset(x, select = c("Year", "count", "norm.count")))
    complete <- expand.grid(year = range[1]:range[2])
    new.vars <- c(paste(var.name, ".count", sep = ""), paste(var.name, ".norm", sep = ""))
    names(x) <- c("year", new.vars)
    x <- merge(complete, x, by = "year", all.x = TRUE)
    x[is.na(x[, paste(var.name, ".norm", sep = "")]), new.vars] <- c(0, 0)
  }
  
  row.names(x) <- 1:nrow(x)
  return(x)
}
