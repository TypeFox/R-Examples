plot.rtconnect <-
function(x, type="daily", graph.type="b",
         main=NULL, xlab=NULL, ylab=NULL, ...) {
  rtc <- x

  if (type == "daily") {
    func.daily <- function(date) {
      sum(subset(rtc, date=date)$Units)
    }

    start <- min(rtc$Begin.Date)
    end <- max(rtc$End.Date)
    days <- as.Date(start:end, origin="1970-01-01")
    daily.units <- sapply(days, func.daily)

    x <- days
    y <- cumsum(daily.units)
    t.main <- sprintf("# of total units from %s to %s",
                      as.character(start), as.character(end))
    t.xlab <- "Day"
    t.ylab <- "# of units"
  } else if (type == "weekly") {
    func.weekly <- function(week) {
      sum(subset(rtc, date.gte=week.to.date(week),
                      date.lte=week.to.date(week)+6)$Units)
    }

    start <- min(rtc$Begin.Date)
    end <- max(rtc$End.Date)
    days <- as.Date(start:end, origin="1970-01-01")
    weeks <- unique(sapply(days, date.to.week))
    weekly.units <- sapply(weeks, func.weekly)

    x <- sapply(weeks, week.to.date)
    y <- cumsum(weekly.units)
    t.main <- sprintf("# of total units from %s to %s",
                      as.character(start), as.character(end))
    t.xlab <- "Week"
    t.ylab <- "# of units"
  }

  if (missing(main)) main <- t.main
  if (missing(xlab)) xlab <- t.xlab
  if (missing(ylab)) ylab <- t.ylab

  par(xaxt="n")
  plot(
    x,
    y,
    type=graph.type,
    main=main,
    xlab=xlab,
    ylab=ylab
  )
  par(xaxt="s")

  if (type == "daily") {
    axis.Date(1, at=seq(start, end, "days"), format="%Y-%m-%d")
  } else if (type == "weekly") {
    axis.Date(1, at=seq(start, end, "week"), format="%YW%U")
  }
}
