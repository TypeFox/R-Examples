origin.year <- 1970
origin.year.POSIXlt <- 1900
class.list <- c("PCICt")

setOldClass("PCICt")

## TODO:
## - S4 class to avoid stripping of attributes?

PCICt.get.months <- function(cal) {
  m.365 <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  m.360 <- c(30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30)
  switch(cal, "365"=m.365, "360"=m.360)
}

dpy.for.cal <- function(cal) {
  switch(cal, "365"=365, "360"=360)
}

clean.cal <- function(cal) {
  cal.list <- c("365_day", "365", "noleap", "360_day", "360", "gregorian",  "standard", "proleptic_gregorian")
  cal.map <- c( "365",     "365", "365",    "360",     "360", "gregorian", "gregorian", "proleptic_gregorian")
  if(!cal %in% cal.list) stop(paste("Calendar type not one of", paste(cal.list, sep=", ")))
  return(cal.map[cal.list %in% cal])
}

.PCICt <- function(x, cal) {
  if(missing(cal)) stop("Can't create a PCICt with no calendar type")
  cal.cleaned <- clean.cal(cal)
  structure(x, cal=cal.cleaned, months=PCICt.get.months(cal.cleaned), class=class.list, dpy=dpy.for.cal(cal.cleaned), tzone="GMT", units="secs")
}

range.PCICt <- function(..., na.rm=FALSE) {
  args <- list(...)
  stopifnot(length(unique(lapply(args, function(x) { attr(x, "cal") }))) == 1)
  args.flat <- unlist(args)
  ret <- c(min(args.flat, na.rm=na.rm), max(args.flat, na.rm=na.rm))
  ret <- copy.atts.PCICt(args[[1]], ret)
  class(ret) <- class.list
  return(ret)
}

c.PCICt <- function(..., recursive=FALSE) {
  ##stopifnot(length(unique(lapply(..., function(x) { attr(x, "cal") }))) == 1)
  cal <- attr(..1, "cal")
  .PCICt(c(unlist(lapply(list(...), unclass))), cal)
}

## Use this to drop the 'units' attribute and unclass the object...
coerceTimeUnit <- function(x) {
  as.vector(switch(attr(x,"units"),
                   secs = x, mins = 60*x, hours = 60*60*x,
                   days = 60*60*24*x, weeks = 60*60*24*7*x))
}

`+.PCICt` <- function(e1, e2) {
  if (nargs() == 1) return(e1)
  ## only valid if one of e1 and e2 is a scalar/difftime
  if(inherits(e1, "PCICt") && inherits(e2, "PCICt"))
    stop("binary '+' is not defined for \"PCICt\" objects")
  if (inherits(e1, "difftime")) e1 <- coerceTimeUnit(e1)
  if (inherits(e2, "difftime")) e2 <- coerceTimeUnit(e2)
  .PCICt(unclass(e1) + unclass(e2), cal=attr(e1, "cal"))
}

`-.PCICt` <- function(e1, e2) {
  ## need to drop "units" attribute here
  if(!inherits(e1, "PCICt"))
    stop("Can only subtract from PCICt objects")
  if (nargs() == 1) stop("unary '-' is not defined for \"PCICt\" objects")
  if(inherits(e2, "PCICt")) {
    stopifnot(attr(e1, "cal") == attr(e2, "cal"))
    return(as.difftime(unclass(e1) - unclass(e2), units="secs"))
  }
  if (inherits(e2, "difftime")) e2 <- coerceTimeUnit(e2)
  if(!is.null(attr(e2, "class")))
    stop("can only subtract numbers from PCICt objects")
  .PCICt(unclass(e1) - e2, cal=attr(e1, "cal"))
}

Ops.PCICt <- function(e1, e2) {
  if (nargs() == 1)
    stop(gettextf("unary '%s' not defined for \"PCICt\" objects",
                  .Generic), domain = NA)

  PCICt.object <- NULL
  if(inherits(e1, "PCICt"))
    PCICt.object <- e1
  else if(inherits(e2, "PCICt"))
    PCICt.object <- e2
  else
    stop("Can't use PCICt operators on non-PCICt objects")


  boolean <- switch(.Generic, "<" = , ">" = , "==" = ,
                    "!=" = , "<=" = , ">=" = TRUE, FALSE)
  if (!boolean)
    stop(gettextf("'%s' not defined for \"PCICt\" objects", .Generic),
         domain = NA)
  if(inherits(e1, "POSIXlt") || is.character(e1)) e1 <- as.PCICt(e1, cal=attr(PCICt.object, "cal"))
  if(inherits(e2, "POSIXlt") || is.character(e1)) e2 <- as.PCICt(e2, cal=attr(PCICt.object, "cal"))
  stopifnot(attr(e1, "cal") == attr(e2, "cal"))
  NextMethod(.Generic)
}

rep.PCICt <- function(x, ...) {
  y <- rep(unclass(x), ...)
  .PCICt(y, cal=attr(x, "cal"))
}

mean.PCICt <- function(x, ...) {
  .PCICt(mean(unclass(x), ...), attr(x, "cal"))
}

min.PCICt <- function(x, ...) {
  res <- min(unclass(x), ...)
  return(copy.atts.PCICt(x, res))
}

max.PCICt <- function(x, ...) {
  res <- max(unclass(x), ...)
  return(copy.atts.PCICt(x, res))
}

seq.PCICt <- function(from, to, by, length.out = NULL, along.with = NULL, ...) {
  if (missing(from))
    stop("'from' must be specified")
  if (!inherits(from, "PCICt"))
    stop("'from' must be a PCICt object")
  if (length(from) != 1L)
    stop("'from' must be of length 1")
  if (!missing(to)) {
    stopifnot(attr(from, "cal") == attr(to, "cal"))
    if (!inherits(to, "PCICt"))
      stop("'to' must be a PCICt object")
    if (length(to) != 1)
      stop("'to' must be of length 1")
    if(to < from)
      stop("'to' must be less than 'from'")
  }
  if (!missing(along.with)) {
    length.out <- length(along.with)
  }
  else if (!is.null(length.out)) {
    if (length(length.out) != 1L)
      stop("'length.out' must be of length 1")
    length.out <- ceiling(length.out)
  }
  status <- c(!missing(to), !missing(by), !is.null(length.out))
  if (sum(status) != 2L)
    stop("exactly two of 'to', 'by' and 'length.out' / 'along.with' must be specified")
  if (missing(by)) {
    from <- unclass(from)
    to <- unclass(to)
    res <- seq.int(from, to, length.out = length.out)
    return(.PCICt(res, attr(from, "cal")))
  }
  if (length(by) != 1L)
    stop("'by' must be of length 1")
  valid <- 0L
  if (inherits(by, "difftime")) {
    by <- switch(attr(by, "units"), secs = 1, mins = 60,
                 hours = 3600, days = 86400, weeks = 7 * 86400) *
                   unclass(by)
  } else if (is.character(by)) {
    by2 <- strsplit(by, " ", fixed = TRUE)[[1L]]
    if (length(by2) > 2L || length(by2) < 1L)
      stop("invalid 'by' string")
    valid <- pmatch(by2[length(by2)], c("secs", "mins", "hours",
                                        "days", "weeks", "months", "years", "DSTdays"))
    if (is.na(valid))
      stop("invalid string for 'by'")
    if (valid <= 5L) {
      by <- c(1, 60, 3600, 86400, 7 * 86400)[valid]
      if (length(by2) == 2L)
        by <- by * as.integer(by2[1L])
    }
    else by <- if (length(by2) == 2L)
      as.integer(by2[1L])
    else 1
  }
  else if (!is.numeric(by))
    stop("invalid mode for 'by'")
  if (is.na(by))
    stop("'by' is NA")
  if (valid <= 5L) {
    from <- unclass(from)
    if (!is.null(length.out))
      res <- seq.int(from, by = by, length.out = length.out)
    else {
      to0 <- unclass(to)
      res <- seq.int(0, to0 - from, by) + from
    }
    return(.PCICt(res, attr(from, "cal")))
  } else {
    r1 <- as.POSIXlt(from)
    if (valid == 7L) {
      if (missing(to)) {
        yr <- seq.int(r1$year, by = by, length.out = length.out)
      } else {
        to0 <- as.POSIXlt(to)
        yr <- seq.int(r1$year, to0$year, by)
      }
      r1$year <- yr
    } else if (valid == 6L) {
      if (missing(to)) {
        mon <- seq.int(r1$mon, by = by, length.out = length.out)
      } else {
        to0 <- as.POSIXlt(to)
        mon <- seq.int(r1$mon, 12 * (to0$year - r1$year) +
                       to0$mon, by)
      }
      r1$mon <- mon
    } else if (valid == 8L) {
      if (!missing(to)) {
        length.out <- 2L + floor((unclass(to) -
                                  unclass(from))/86400)
      }
      r1$mday <- seq.int(r1$mday, by = by, length.out = length.out)
    }
    r1$isdst <- -1L
    res <- as.PCICt(r1, attr(from, "cal"))
    if (!missing(to)) {
      res <- if (by > 0)
        res[res <= to]
      else res[res >= to]
    }
    res
  }
}

trunc.PCICt <- function(x, units = c("secs", "mins", "hours", "days"), ...) {
      units <- match.arg(units)
      val <- unclass(x)
      round.to <- switch(units, secs = 1, mins = 60, hours = 3600, days = 86400)
      val <- floor(val / round.to) * round.to
      class(val) <- class(x)
      return(copy.atts.PCICt(x, val))
}

round.PCICt <- function (x, digits = c("secs", "mins", "hours", "days")) {
  if (is.numeric(digits) && digits == 0)
    digits <- "secs"
  digits <- match.arg(digits)
  x <- x + switch(digits, secs = 0.5, mins = 30, hours = 1800,
                  days = 43200)
  trunc(x, units = digits)
}

copy.atts.PCICt <- function(from, to) {
  return(structure(to, cal=attr(from, "cal"), months=attr(from, "months"), class=class(from), dpy=attr(from, "dpy"), tzone=attr(from, "tzone"), units=attr(from, "units")))
}

`[.PCICt` <- function(x, ...) {
  val <- NextMethod("[")
  val <- copy.atts.PCICt(x, val)
  class(val) <- class(x)
  val
}

`[<-.PCICt` <- function (x, ..., value) {
  if (!as.logical(length(value)))
    return(x)
  stopifnot(class(value) == class(x) & attr(x, "cal") == attr(value, "cal"))
  cl <- oldClass(x)
  x <- NextMethod("[<-")
  x <- copy.atts.PCICt(value, x)
  class(x) <- cl
  x
}

as.PCICt <- function(x, cal, ...) {
  if(missing(cal)) stop("Can't create a PCICt with no calendar type")
  UseMethod("as.PCICt")
}

as.character.PCICt <- function(x, ...) {
  format.PCICt(x, ...)
}

unique.PCICt <- function(x, incomparables = FALSE, fromLast = FALSE, ...) {
  if (!inherits(x, "PCICt"))
    stop("wrong class")
  z <- unique(unclass(x), incomparables, fromLast, ...)
  return(copy.atts.PCICt(x, z))
}

summary.PCICt <- function (object, digits = 15, ...) {
  x <- summary.default(unclass(object), digits = digits, ...)
  if (m <- match("NA's", names(x), 0)) {
    NAs <- as.integer(x[m])
    x <- x[-m]
    attr(x, "NAs") <- NAs
  }
  x <- copy.atts.PCICt(object, x)
  class(x) <- c("summaryDefault", "table", oldClass(object))
  x
}

format.PCICt <- function(x, format="", tz="", usetz=FALSE, ...) {
  if (!inherits(x, "PCICt"))
    stop("wrong class")
  
  if(!is.null(attr(x, "dpy")) && attr(x, "dpy") == 360) {
    structure(format.POSIXlt.360(as.POSIXlt(x, tz), format,
                             ...), names = names(x))
  } else {
    structure(format.POSIXlt(as.POSIXlt(x, tz), format, usetz,
                             ...), names = names(x))
  }
}

print.PCICt <- function (x, ...) {
  max.print <- getOption("max.print", 9999L)
  if (max.print < length(x)) {
    print(as.character(x[1:max.print]), ...)
    cat(" [ reached getOption(\"max.print\") -- omitted",
        length(x) - max.print, "entries ]\n")
  }
  else print(as.character(x), ...)
  invisible(x)
}


strptime.360 <- function(x, format) {
  .Call("do_strptime_360", x, format)
}

format.POSIXlt.360 <- function(x, format="") {
  if (!inherits(x, "POSIXlt"))
    stop("wrong class")
  if (format == "") {
    times <- unlist(unclass(x)[1L:3L])
    secs <- x$sec
    secs <- secs[!is.na(secs)]
    np <- getOption("digits.secs")
    if (is.null(np))
      np <- 0L
    else np <- min(6L, np)
    if (np >= 1L)
      for (i in seq_len(np) - 1L)
        if (all(abs(secs - round(secs, i)) < 1e-06)) {
        np <- i
        break
      }
    format <- if (all(times[!is.na(times)] == 0))
      "%Y-%m-%d"
    else if (np == 0L)
      "%Y-%m-%d %H:%M:%S"
    else paste("%Y-%m-%d %H:%M:%OS", np, sep="")
  }
  y <- .Call("do_formatPOSIXlt_360", x, format)
  names(y) <- names(x$year)
  y
  
}

as.POSIXct.POSIXlt.360 <- function(x) {
  .Call("do_asPOSIXct_360", x, format)
}

as.POSIXlt.POSIXct.360 <- function(x) {
  .Call("do_asPOSIXlt_360", x, format)
}

as.PCICt.default <- function(x, cal, format, ...) {
  tz <- "GMT"
  cal.cleaned <- clean.cal(cal)
  if (inherits(x, "PCICt"))
    return(x)
  if (is.character(x) || is.factor(x)) {
    x <- as.character(x)
    if(cal.cleaned == "360") {
      if (!missing(format)) {
        res <- strptime.360(x, format)
        return(as.PCICt(res, cal, ...))
      }
      x <- unclass(x)
      xx <- x[!is.na(x)]
      if (!length(xx))
        res <- strptime.360(x, "%Y/%m/%d")
      else if (all(!is.na(strptime.360(xx, f <- "%Y-%m-%d %H:%M:%OS"))) ||
               all(!is.na(strptime.360(xx, f <- "%Y/%m/%d %H:%M:%OS"))) ||
               all(!is.na(strptime.360(xx, f <- "%Y-%m-%d %H:%M"))) ||
               all(!is.na(strptime.360(xx, f <- "%Y/%m/%d %H:%M"))) ||
               all(!is.na(strptime.360(xx, f <- "%Y-%m-%d"))) ||
               all(!is.na(strptime.360(xx, f <- "%Y/%m/%d"))))
        res <- strptime.360(x, f)
      if(missing(res)) stop("character string is not in a standard unambiguous format")
      return(as.PCICt(res, cal, ...))
    } else {
      return(as.PCICt(as.POSIXlt(x, tz, format, ...), cal, ...))
    }
  }
  if (is.logical(x) && all(is.na(x)))
    return(.PCICt(as.numeric(x), cal))
  stop(gettextf("do not know how to convert '%s' to class \"PCICt\"", deparse(substitute(x))))
}

as.PCICt.numeric <- function(x, cal, origin, ...) {
  if (missing(origin))
    stop("'origin' must be supplied")

  if(inherits(origin, "PCICt") && attr(origin, "cal") == cal)
    return(origin + x)
  else
    return(as.PCICt(origin, cal) + x)
}

as.PCICt.POSIXlt <- function(x, cal, ...) {
  proleptic.correction <- 0
  seconds.per.day <- 86400
  tz <- "GMT"
  cal.cleaned <- clean.cal(cal)
  year.length <- dpy.for.cal(cal.cleaned)

  if(is.null(year.length)) {
    d <- as.POSIXct(x, tz="GMT")
    class(d) <- NULL
    return(.PCICt(d, "proleptic_gregorian"))
  } else {
    months <- PCICt.get.months(cal.cleaned)
    months.off <- cumsum(c(0, months[1:(length(months) - 1)]))
    seconds.per.hour <- 3600
    return(.PCICt((x$year + origin.year.POSIXlt - origin.year + floor(x$mon / 12)) * year.length * seconds.per.day +
                  months.off[(x$mon %% 12) + 1] * seconds.per.day + (x$mday - 1) * seconds.per.day + x$hour * seconds.per.hour + x$min * 60 + x$sec, cal=cal))
  }
}

as.PCICt.POSIXct <- function(x, cal, ...) {
  cal.cleaned <- clean.cal(cal)
  if(cal.cleaned == "360") {
    as.PCICt.POSIXlt(as.POSIXlt.POSIXct.360(x), cal, ...)
  } else {
    as.PCICt.POSIXlt(as.POSIXlt(x), cal, ...)
  }
}

## FIXME: Better NA handling
as.POSIXlt.PCICt <- function(x, tz="", ...) {
  seconds.per.day <- 86400
  seconds.per.hour <- 3600

  tzone <- attr(x, "tzone")
  if (length(tz) == 0 && !is.null(tzone))
    tz <- tzone[1L]

  if(is.null(attr(x, "months"))) {
    class(x) <- c("POSIXct", "POSIXt")
    return(as.POSIXlt(x))
  } else {
    months <- attr(x, "months")
    months.off <- cumsum(c(0, months[1:(length(months) - 1)]))
    months.idx <- unlist(lapply(1:12, function(x) { rep(x, months[x]) } ))

    days.per.year <- attr(x, "dpy")
    remainder <- as.numeric(x) %% (days.per.year * seconds.per.day)
    remainder[remainder < 0] <- days.per.year * seconds.per.day - remainder[remainder < 0]

    year <- floor(as.numeric(x) / (days.per.year * seconds.per.day)) + origin.year
    yday <- floor(remainder / seconds.per.day) + 1
    month <- months.idx[yday]
    day <- yday - months.off[month]

    ## Need to compute wday
    wday <- (as.numeric(x) / 86400) %% 7
    hms.remainder <- remainder %% seconds.per.day
    hour <- floor(hms.remainder / seconds.per.hour)
    minute <- floor((hms.remainder %% seconds.per.hour) / 60)
    second <- hms.remainder %% 60
    return(.POSIXlt(list(sec=second, min=minute, hour=hour, mday=day, mon=month - 1, year=year - origin.year.POSIXlt, wday=wday, yday=yday - 1, isdst=0), tz))
  }
}

as.POSIXct.PCICt <- function(x, tz="", ...) {
  
  if(attr(x, "cal") == "360") {
    warning("360-day PCICt objects can't be properly represented by a POSIXct object")
  }
  return(as.POSIXct(as.POSIXlt(x, tz)))
}

cut.PCICt <- function (x, breaks, labels = NULL, start.on.monday = TRUE, right = FALSE, ...) {
  if(!inherits(x, "PCICt")) stop("'x' must be a PCICt object")
  cal <- attr(x, "cal")
  
  if (inherits(breaks, "PCICt") || (is.numeric(breaks) && length(breaks) == 1L)) {
    ## Dates are already PCICt or specified number of breaks; don't need to do anything
  } else if(is.character(breaks) && length(breaks) == 1L) {
    ## Breaks are characters; need to do something.
    by2 <- strsplit(breaks, " ", fixed=TRUE)[[1L]]
    if(length(by2) > 2L || length(by2) < 1L)
      stop("invalid specification of 'breaks'")
    valid <- pmatch(by2[length(by2)],
                    c("secs", "mins", "hours", "days", "weeks",
                      "months", "years", "DSTdays", "quarters"))
    if(is.na(valid)) stop("invalid specification of 'breaks'")
    start <- as.POSIXlt(min(x, na.rm=TRUE))
    incr <- 1
    if(valid > 1L) { start$sec <- 0L; incr <- 60 }
    if(valid > 2L) { start$min <- 0L; incr <- 3600 }
    ## start of day need not be on the same DST, PR#14208
    if(valid > 3L) { start$hour <- 0L; start$isdst <- -1L; incr <- 86400 }
    if(valid == 5L) {               # weeks
      start$mday <- start$mday - start$wday
      if(start.on.monday)
        start$mday <- start$mday + ifelse(start$wday > 0L, 1L, -6L)
      incr <- 7*86400
    }
    if(valid == 8L) incr <- 25*3600 # DSTdays
    if(valid == 6L) {               # months
      start$mday <- 1L
      maxx <- max(x, na.rm = TRUE)
      step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
      end <- as.POSIXlt(maxx + (ifelse(cal == "360", 30, 31) * step * 86400))
      end$mday <- 1L
      end$isdst <- -1L
      breaks <- seq(as.PCICt(start, cal), as.PCICt(end, cal), breaks)
    } else if(valid == 7L) {        # years
      start$mon <- 0L
      start$mday <- 1L
      maxx <- max(x, na.rm = TRUE)
      step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
      end <- as.POSIXlt(maxx + (ceiling(get.avg.dpy(x)) * step* 86400))
      end$mon <- 0L
      end$mday <- 1L
      end$isdst <- -1L
      breaks <- seq(as.PCICt(start, cal), as.PCICt(end, cal), breaks)
    } else if(valid == 9L) {        # quarters
      qtr <- rep(c(0L, 3L, 6L, 9L), each = 3L)
      start$mon <- qtr[start$mon + 1L]
      start$mday <- 1L
      maxx <- max(x, na.rm = TRUE)
      step <- ifelse(length(by2) == 2L, as.integer(by2[1L]), 1L)
      end <- as.POSIXlt(maxx + (floor(get.avg.dpy(x) / 4) * step * 86400))
      end$mon <- qtr[end$mon + 1L]
      end$mday <- 1L
      end$isdst <- -1L
      breaks <- seq(as.PCICt(start, cal), as.PCICt(end, cal), paste(step * 3, "months"))
      ## 90-93 days ahead could give an empty level, so
      lb <- length(breaks)
      if(maxx < breaks[lb-1]) breaks <- breaks[-lb]
    } else {                        # weeks or shorter
      if (length(by2) == 2L) incr <- incr * as.integer(by2[1L])
      maxx <- max(x, na.rm = TRUE)
      breaks <- seq(as.PCICt(start, cal), maxx + incr, breaks)
      breaks <- breaks[seq_len(1+max(which(breaks <= maxx)))]
    }
  } else stop("invalid specification of 'breaks'")
  res <- cut(unclass(x), unclass(breaks), labels = labels,
             right = right, ...)
  if(is.null(labels)) {
    levels(res) <-
      as.character(if (is.numeric(breaks)) x[!duplicated(res)]
      else breaks[-length(breaks)])
  }
  res
}

diff.PCICt <- function (x, lag = 1L, differences = 1L, ...) {
  class(x) <- c("POSIXct", "POSIXt")
  diff(x, lag, differences, ...)
}

is.numeric.PCICt <- function(x) FALSE

julian.PCICt <- function (x, origin=NULL, ...) {
  if(is.null(origin))
    origin <- "1970-01-01"
  else
    stopifnot(attr(x, "cal") == attr(origin, "cal"))

  origin <- as.PCICt(origin, cal=attr(x, "cal"))
  class(x) <- class(origin) <- c("POSIXct", "POSIXt")
  if (length(origin) != 1L)
    stop("'origin' must be of length one")

  res <- difftime(x, origin, units = "days")
  structure(res, origin = origin)
}

get.sec.incr <- function(x, secs, incr=1, mul=1.1) {
  if(length(secs) == 0 || mul * (incr * secs[1]) > x)
    incr
  else
    get.sec.incr(x, secs[-1], incr * secs[1], mul)
}

Axis.PCICt <- function(x = NULL, at = NULL, ..., side, labels = TRUE) {
  axis.PCICt(side = side, x = x, at = at, labels = labels, ...)
}

get.avg.dpy <- function(x) {
  ifelse(is.null(attr(x, "dpy")), 365.25, attr(x, "dpy"))
}

axis.PCICt <- function(side, x, at, format, labels = TRUE, ...) {
  mat <- missing(at) || is.null(at)
  mft <- missing(format) || is.null(format)
  if (!mat)
    x <- at

  range <- par("usr")[if (side%%2) 1L:2L else 3L:4L]

  d <- range[2L] - range[1L]
  z <- c(as.PCICt(range, cal=attr(x, "cal"), origin="1970-01-01"), x[is.finite(x)])

  sc <- get.sec.incr(d, c(60, 60, 24, 7))
  if(mft && !is.na(sc))
    format <- switch(as.character(sc), "1"="%S", "60"="%M:%S", "3600"="%H:%M", "86400"="%a %H:%M", "604800"="%a")

  if (d < 60 * 60 * 24 * 50) {
    zz <- pretty(unclass(z)/sc)
    z <- .PCICt(zz * sc, cal=attr(x, "cal"))
    if (!is.na(sc) && sc == 60 * 60 * 24)
      z <- round(z, "days")
    if (mft)
      format <- "%b %d"
  } else if (d < 1.1 * 60 * 60 * 24 * get.avg.dpy(x)) {
    zz <- as.POSIXlt(z)
    zz$mday <- zz$wday <- zz$yday <- 1
    zz$isdst <- -1
    zz$hour <- zz$min <- zz$sec <- 0
    zz$mon <- pretty(zz$mon)
    m <- length(zz$mon)
    M <- 2 * m
    m <- rep.int(zz$year[1L], m)
    zz$year <- c(m, m + 1)
    zz <- lapply(zz, function(x) rep(x, length.out = M))
    z <- as.PCICt(zz, attr(x, "cal"))
    if (mft)
      format <- "%b"
  } else {
    zz <- as.POSIXlt(z)
    zz$mday <- zz$wday <- zz$yday <- 1
    zz$isdst <- -1
    zz$mon <- zz$hour <- zz$min <- zz$sec <- 0
    zz$year <- pretty(zz$year)
    M <- length(zz$year)
    zz <- lapply(zz, function(x) rep(x, length.out = M))
    z <- as.PCICt(.POSIXlt(zz), attr(x, "cal"))
    if (mft)
      format <- "%Y"
  }
  if (!mat)
    z <- x[is.finite(x)]

  keep <- z >= range[1L] & z <= range[2L]
  z <- z[keep]
  if (!is.logical(labels))
    labels <- labels[keep]
  else if (identical(labels, TRUE))
    labels <- format(z, format = format)
  else if (identical(labels, FALSE))
    labels <- rep("", length(z))
  
  axis(side, at = unclass(z), labels = labels, ...)
}

pretty.PCICt <- function(x, n = 5, min.n = n %/% 2, sep = " ", ...) {
  zz <- range(x, na.rm = TRUE)
  xspan <- as.numeric(diff(zz), units = "secs")
  if (diff(as.numeric(zz)) == 0) # one value only
    zz <- zz + c(0,60)
  ## specify the set of pretty timesteps
  MIN <- 60
  HOUR <- MIN * 60
  DAY <- HOUR * 24
  YEAR <- DAY * get.avg.dpy(x)
  MONTH <- YEAR / 12
  steps <-
    list("1 sec" = list(1, format = "%S", start = "mins"),
         "2 secs" = list(2),
         "5 secs" = list(5),
         "10 secs" = list(10),
         "15 secs" = list(15),
         "30 secs" = list(30, format = "%H:%M:%S"),
         "1 min" = list(1*MIN, format = "%H:%M"),
         "2 mins" = list(2*MIN, start = "hours"),
         "5 mins" = list(5*MIN),
         "10 mins" = list(10*MIN),
         "15 mins" = list(15*MIN),
         "30 mins" = list(30*MIN),
         ## "1 hour" = list(1*HOUR),
         "1 hour" = list(1*HOUR, format = if (xspan <= DAY) "%H:%M" else paste("%b %d", "%H:%M", sep = sep)),
         "3 hours" = list(3*HOUR, start = "days"),
         "6 hours" = list(6*HOUR, format = paste("%b %d", "%H:%M", sep = sep)),
         "12 hours" = list(12*HOUR),
         "1 DSTday" = list(1*DAY, format = paste("%b", "%d", sep = sep)),
         "2 DSTdays" = list(2*DAY),
         "1 week" = list(7*DAY, start = "weeks"),
         "halfmonth" = list(MONTH/2, start = "months"),
         ## "1 month" = list(1*MONTH, format = "%b"),
         "1 month" = list(1*MONTH, format = if (xspan < YEAR) "%b" else paste("%b", "%Y", sep = sep)),
         "3 months" = list(3*MONTH, start = "years"),
         "6 months" = list(6*MONTH, format = "%Y-%m"),
         "1 year" = list(1*YEAR, format = "%Y"),
         "2 years" = list(2*YEAR, start = "decades"),
         "5 years" = list(5*YEAR),
         "10 years" = list(10*YEAR),
         "20 years" = list(20*YEAR, start = "centuries"),
         "50 years" = list(50*YEAR),
         "100 years" = list(100*YEAR),
         "200 years" = list(200*YEAR),
         "500 years" = list(500*YEAR),
         "1000 years" = list(1000*YEAR))
  ## carry forward 'format' and 'start' to following steps
  for (i in seq_along(steps)) {
    if (is.null(steps[[i]]$format))
      steps[[i]]$format <- steps[[i-1]]$format
    if (is.null(steps[[i]]$start))
      steps[[i]]$start <- steps[[i-1]]$start
    steps[[i]]$spec <- names(steps)[i]
  }
  ## crudely work out number of steps in the given interval
  nsteps <- sapply(steps, function(s) {
    xspan / s[[1]]
  })
  init.i <- which.min(abs(nsteps - n))
  ## calculate actual number of ticks in the given interval
  calcSteps <- function(s) {
    startTime <- trunc(min(zz), units = s$start)
    if (identical(s$spec, "halfmonth")) {
      at <- seq(startTime, max(zz), by = "months")
      at2 <- as.POSIXlt(at)
      at2$mday <- 15L
      at3 <- sort(c(at, as.PCICt(at2)))
      at <- copy.atts.PCICt(at, at3)
    } else {
      at <- seq(startTime, max(zz), by = s$spec)
    }
    at <- at[(min(zz) <= at) & (at <= max(zz))]
    at
  }
  init.at <- calcSteps(steps[[init.i]])
  init.n <- length(init.at) - 1L
  ## bump it up if below acceptable threshold
  while (init.n < min.n) {
    init.i <- init.i - 1L
    if (init.i == 0) stop("range too small for min.n")
    init.at <- calcSteps(steps[[init.i]])
    init.n <- length(init.at) - 1L
  }
  makeOutput <- function(at, s) {
    flabels <- format(at, s$format)
    ans <- as.PCICt(at, cal=attr(x, "cal"))
    attr(ans, "labels") <- flabels
    ans
  }
  if (init.n == n) ## perfect
    return(makeOutput(init.at, steps[[init.i]]))
  if (init.n > n) {
    ## too many ticks
    new.i <- init.i + 1L
    new.i <- min(new.i, length(steps))
  } else {
    ## too few ticks
    new.i <- init.i - 1L
    new.i <- max(new.i, 1L)
  }
  new.at <- calcSteps(steps[[new.i]])
  new.n <- length(new.at) - 1L
  ## work out whether new.at or init.at is better
  if (new.n < min.n)
    new.n <- -Inf
  if (abs(new.n - n) < abs(init.n - n))
    makeOutput(new.at, steps[[new.i]])
  else
    makeOutput(init.at, steps[[init.i]])
}
