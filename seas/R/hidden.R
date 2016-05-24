".seasxlab" <- 
function(width, start.day) {
  str <- if (is.numeric(width))
    gettextf("%s-day group", round(width, 1))
  else if (any(grep("mon", width)))
    gettext("Monthly")
  else if (any(grep("zod", width)))
    gettext("Zodiac")
  else if (width == "DJF")
    gettext("Quarterly (seasonal)")
  else if (width == "JFM")
    gettext("Quarterly (annual)")
  else if (width == "JF")
    gettext("Two months")
  else as.character(width[1])
  if (!is.null(start.day) && start.day != 1) {
    if (inherits(start.day, c("POSIXct", "Date"))) {
      mday <- as.integer(format(start.day, "%d"))
      fmt <- if (mday == 1)
        getOption("seas.label")$month
      else
        getOption("seas.label")$monthday
      st <- format(start.day, fmt)
    } else {
      mday <- 2
      st <- gettextf("day %i of the year", start.day)
    }
    str <- gettextf(ngettext(mday,
                             "%s, starting in %s",
                             "%s, starting on %s"), str, st)
  }
  return(str)
}

".seasylab" <-
function(var, long.name = NULL, units = NULL) {
  ylab <- if (is.null(long.name)) var else long.name
  if (!is.null(units))
    ylab <- sprintf(getOption("seas.label")$fmt, ylab, units)
  ylab
}

".seastitle" <- 
function(main = NULL, id = NULL, name = NULL,
         orig = NULL, fun = NULL, range = NA) {
  op <- getOption("seas.main")
  if (is.null(main)) {
    if (is.null(name) && !is.null(id) && op$show.id)
      name <- getstnname(id)
    if (!is.null(name) && !is.null(id) && op$show.id)
      main <- paste(name, id)
    else if (!is.null(name) && (!op$show.id || is.null(id)))
      main <- name
    else if (!is.null(id) && op$show.id)
      main <- id
    else if (!is.null(orig))
      main <- orig[[1]]
    else
      main <- ""
  }
  if (!is.null(fun) && op$show.fun)
    main <- paste(main, fun, sep=ifelse(main == "", "", " "))
  if (length(range) == 2 && all(is.finite(range))) {
    range <- ifelse(range[1] == range[2],
                    range[1],
                    paste(range[1], range[2], sep=op$rngsep))
    main <- sprintf(op$fmt, main, range)
  }
  main
}

".seasmonthgrid" <-
function(width, days, start = 1, rep = 0, start.day = 1, month.label) {
  op <- getOption("seas.month.grid")
  if (is.numeric(width)) {
    num.bin <- length(days)
    num <- num.bin + rep
    year.length <- round(sum(days))
    month.days <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    if (year.length != 366) {
      if (year.length == 365)
        month.days[2] <- 28
      else
        month.days <- rep(year.length / 12, 12)
    }
    if (inherits(start.day, c("POSIXct", "Date"))) {
      start.std <- as.Date(format(start.day, "2000-%m-%d"))
      start.yday <- as.integer(format(start.std, "%j"))
    } else {
      start.yday <- start.day
    }
    if (start > 1)
      start.yday <- round(sum(start.yday-1,
                              days[1:(start - 1)])) %% year.length + 1
    start.date <- as.Date(sprintf("2000-%03i", as.integer(start.yday)), "%Y-%j")
    sel.bin <- seq(start - 1, start + num - 2) %% num.bin + 1
    bdf <- data.frame(pos=1:num, bin=sel.bin, width=days[sel.bin])
    bdf$yday <- cumsum(bdf$width) - bdf$width + start.yday
    yday2pos <- function(yday) {  # y=mx+b, using intervals
      p <- findInterval(yday, bdf$yday)
      if (p < 1) p <- 1
      m <- 1 / bdf$width[p]
      b <- bdf$pos[p] - m * bdf$yday[p]
      return(m * yday + b - 0.5)
    }
    ym <- as.numeric(format(start.date, "%m"))
    yd <- ifelse(ym > 1, sum(month.days[1:(ym - 1)]), 1)
    pos <- yday2pos(yd)
    ln <- 1
    while (pos[ln] < num) {
      ym[ln + 1] <- ym[ln] %% 12 + 1
      yd <- yd + month.days[ym[ln + 1]]
      pos[ln + 1] <- yday2pos(yd)
      ln <- length(ym)
    }
    fit <- pos < num + 1.4
    pos <- pos[fit]
    ym <- ym[fit]
    month <- months(as.Date(paste(2001, ym[-length(ym)],
                                  15, sep="-")), op$abb)
    if (is.numeric(op$len))
      month <- strtrim(month, op$len)
    abline(v=pos, col=op$col, lwd=op$lwd * par("lwd"), lty=op$lty)
    m <- pos[-length(ym)] + diff(pos) / 2
    fit <- m > 0.2
    m <- m[fit]
    month <- month[fit]
    if (missing(month.label))
      month.label <- op$label
    if (month.label) {
      if (op$force)
        mtext(month, at=m, cex=par("cex") * par("cex.axis"))
      else
        axis(3, m, month, tick=FALSE, line=-1)
    }
  }
}
