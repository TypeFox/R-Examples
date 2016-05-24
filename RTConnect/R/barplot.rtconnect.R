barplot.rtconnect <-
function(height, type="daily", main=NULL, xlab=NULL, ylab=NULL,
         plab=FALSE, ...) {
  rtc <- height
  names.arg <- NULL
  legend.text <- NULL
  col <- NULL

  if (type == "daily") {
    daily <- function(date) {
      sum(subset(rtc, date=date)$Units)
    }

    start <- min(rtc$Begin.Date)
    end <- max(rtc$End.Date)
    days <- as.Date(start:end, origin="1970-01-01")

    units <- sapply(days, daily)
    names.arg <- days
    t.main <- sprintf("# of units from %s to %s",
                   as.character(start), as.character(end))
    t.xlab <- "Day"
    t.ylab <- "# of units"
  } else if (type == "daily.version") {
    daily.version <- function(date) {
      version <- function(version) {
        sum(subset(subset(rtc, version=version), date=date)$Units)
      }
      sapply(versions, version)
    }

    start <- min(rtc$Begin.Date)
    end <- max(rtc$End.Date)
    dates <- as.Date(start:end, origin="1970-01-01")
    versions <- sort(unique(rtc$Version))

    units <- sapply(dates, daily.version)
    names.arg=dates
    t.main <- sprintf("# of daily units from %s to %s",
                   as.character(start), as.character(end))
    t.xlab <- "Date"
    t.ylab <-"# of units"
    legend.text <- versions
    col <- rainbow(nrow(units))
  } else if (type == "weekly") {
    weekly <- function(week) {
      sum(subset(rtc, date.gte=week.to.date(week),
                      date.lte=week.to.date(week)+6)$Units)
    }

    start <- min(rtc$Begin.Date)
    end <- max(rtc$End.Date)
    days <- as.Date(start:end, origin="1970-01-01")
    weeks <- unique(sapply(days, date.to.week))

    units <- sapply(weeks, weekly)
    names.arg <- weeks
    t.main <- sprintf("# of weekly units from %s to %s",
                     as.character(weeks[1]), as.character(weeks[length(weeks)]))
    t.xlab <- "Week"
    t.ylab <-"# of units"
  } else if (type == "weekly.version") {
    weekly.version <- function(week) {
      version <- function(version) {
        sum(subset(rtc, version=version,
                        date.gte=week.to.date(week),
                        date.lte=week.to.date(week)+6)$Units)
      }
      sapply(versions, version)
    }

    start <- min(rtc$Begin.Date)
    end <- max(rtc$End.Date)
    days <- as.Date(start:end, origin="1970-01-01")
    weeks <- unique(sapply(days, date.to.week))
    versions <- sort(unique(rtc$Version))

    units <- sapply(weeks, weekly.version)
    names.arg <- weeks
    t.main <- sprintf("# of weekly units from %s to %s",
                   as.character(start), as.character(end))
    t.xlab <- "Week"
    t.ylab <-"# of units"
    legend.text <- versions
    col <- rainbow(nrow(units))
  } else if (type == "version") {
    version <- function(version) {
      sum(subset(rtc, version=version)$Units)
    }

    versions <- sort(unique(rtc$Version))

    units <- sapply(versions, version)
    names.arg <- versions
    t.main <- sprintf("# of units in %s versions", as.character(length(versions)))
    t.xlab <- "Version"
    t.ylab <-"# of units"
  } else if (type == "country") {
    country <- function(country) {
      sum(subset(rtc, country.code=country)$Units)
    }

    countries <- sort(unique(rtc$Country.Code))

    units <- sapply(countries, country)
    names.arg <- countries
    t.main <- sprintf("# of units in %s countries", as.character(length(countries)))
    t.xlab <- "Country"
    t.ylab <-"# of units"
  }

  if (missing(main)) main <- t.main
  if (missing(xlab)) xlab <- t.xlab
  if (missing(ylab)) ylab <- t.ylab

  bp <- barplot(
    units,
    names.arg=names.arg,
    ylim=c(0, 1.1 * max(units)),
    main=main,
    xlab=xlab,
    ylab=ylab,
    legend.text=legend.text,
    col=col
  )

  if (plab) {
    text(bp, units, units, adj=c(.5,-.5))
  }
}
