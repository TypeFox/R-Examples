"write.sdsm" <-
function(dat, var, start, end, file = "") {
  start <- as.Date(paste(start, 1, 1, sep="-"))
  end <- as.Date(paste(end, 12, 31, sep="-"))
  sdat <- data.frame(date=seq(start, end, by="day"), var=NA)
  dat <- merge(sdat, dat[,c("date", var)], by="date", all.x=TRUE)[,var]
  dat <- sprintf("%9.3f", ifelse(is.na(dat), -999, dat))
  write.table(dat, file, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

"read.sdsm" <-
function(file, start = 1961, end = 2000, calendar) {
  if (!missing(calendar) && !is.null(calendar)) {
    nyears <- end - start + 1
    if (calendar == "360_day") {
      years <- rep(start:end, each=360)
      months <- rep(1:12, nyears, each=30)
      days <- rep.int(rep.int(1:30, 12), nyears)
      days[months == 2 & days == 30] <- 29  # somehow pack 30 valid days in Feb
      days[!(years %% 4 == 0 & years %% 100 !=0 | years %% 400 == 0) & days == 29] <- 28
    } else if (calendar %in% c("365_day", "noleap")) {
      years <- rep(start:end, each=365)
      dpm <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      months <- rep(1:12, times=dpm)
      days <- sequence(dpm)
    } else if (calendar %in% c("366_day", "all_leap", "julian")) {
      years <- rep(start:end, each=366)
      dpm <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      months <- rep(1:12, times=dpm)
      days <- sequence(dpm)
      if (calendar=="julian") {
        valid <- years %% 4 != 0 & months == 2 & days != 29
        years <- years[valid]
        months <- months[valid]
        days <- days[valid]
      }
    } else {
      stop("Unrecognized 'calendar' argument")
    }
    dates <- as.Date(paste(years, months, days, sep="-"))
  } else {
    calendar <- NULL
    start <- as.Date(paste(start, 1, 1, sep="-"))
    end <- as.Date(paste(end, 12, 31, sep="-"))
    dates <- seq(start, end, by="day")
  }
  x <- read.table(file, na.strings="-999.000")
  x$date <- dates
  attr(x$date, "calendar") <- calendar
  x
}
