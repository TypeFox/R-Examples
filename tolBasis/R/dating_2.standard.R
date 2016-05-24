
#-----------------------------------------------------------------------------
# Standard Datings

Yearly <- new_Dating("Yearly")
Monthly <- new_Dating("Monthly")
Weekly <- new_Dating("Weekly")
# Daily <- new_Dating("Daily")

.DDunit <- function(dating) UseMethod(".DDunit")

.DDunit.Yearly <- function(dating) { "year" }
.DDunit.Monthly <- function(dating) { "month" }
.DDunit.Weekly <- function(dating) { "week" }
.DDunit.Daily <- function(dating) { "day" }
.DDunit.default <- function(dating) { "" }

#-----------------------------------------------------------------------------

Dbelong.Dating <- function(dte, dating) {
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  stopifnot(inherits(dating, "Dating"))
  stdunit = .DDunit(dating)
  if(nchar(stdunit)) round_date(dte, stdunit)==dte
  else NextMethod("Dbelong", dating)
}

Dseq.Dating <- function(from, to, dating, len) {
  stopifnot(inherits(dating, "Dating"))
  stdunit = .DDunit(dating)
  if(nchar(stdunit)) {
    stopifnot(inherits(from, c("Date", "POSIXt")))
    from <- as.Date(from)
    if(missing(to)) {
      seq(ceiling_date(from, stdunit), length.out=len, by=stdunit)
    } else {
      stopifnot(inherits(to, c("Date", "POSIXt")))
      to <- as.Date(to)
      stopifnot(from <= to)
      seq(ceiling_date(from, stdunit), floor_date(to, stdunit), stdunit)
    }
  } else NextMethod("Dseq", dating)
}

Dsucc.Dating <- function(dte, dating=Daily, num=1) {
  stopifnot(inherits(dating, "Dating"))
  stdunit = .DDunit(dating)
  if(nchar(stdunit)) {
    stopifnot(inherits(dte, c("Date", "POSIXt")))
    dte <- as.Date(dte)
    stopifnot(inherits(num, c("numeric", "integer")))
    stopifnot(is.finite(num))
    if(num>0) {
      floor_date(dte, stdunit) + period(num, units=stdunit)
    } else if (num<0) {
      ceiling_date(dte, stdunit) - period(abs(num), units=stdunit)
    } else {
      round_date(dte, stdunit)
    }
  } else NextMethod("Dsucc", dating)
}

Dfloor.Dating <- function(dte, dating=Daily) {
  stopifnot(inherits(dating, "Dating"))
  stdunit = .DDunit(dating)
  if(nchar(stdunit)) {
    stopifnot(inherits(dte, c("Date", "POSIXt")))
    dte <- as.Date(dte)
    floor_date(dte, stdunit)
  } else NextMethod("Dfloor", dating)
}

Dceiling.Dating <- function(dte, dating=Daily) {
  stopifnot(inherits(dating, "Dating"))
  stdunit = .DDunit(dating)
  if(nchar(stdunit)) {
    stopifnot(inherits(dte, c("Date", "POSIXt")))
    dte <- as.Date(dte)
    ceiling_date(dte, stdunit)
  } else NextMethod("Dceiling", dating)
}

Dround.Dating <- function(dte, dating=Daily) {
  stopifnot(inherits(dating, "Dating"))
  stdunit = .DDunit(dating)
  if(nchar(stdunit)) {
    stopifnot(inherits(dte, c("Date", "POSIXt")))
    dte <- as.Date(dte)
    round_date(dte, stdunit)
  } else NextMethod("Dround", dating)
}

Ddiff.Dating <- function(dte1, dte2, dating=Daily) {
  stopifnot(inherits(dating, "Dating"))
  #if(as.Date(dte1)>as.Date(dte2)) return(Ddiff.Dating(dte2, dte1, dating))
  NextMethod("Ddiff", dating)
}


#-----------------------------------------------------------------------------
# Specific Methods

Ddiff.Yearly <- function(dte1, dte2, dating=NA) {
  stopifnot(inherits(dte1, c("Date", "POSIXt")))
  stopifnot(inherits(dte2, c("Date", "POSIXt")))
  dt1 <- floor_date(dte1, "year")
  dt2 <- floor_date(dte2, "year")
  year(dt2)-year(dt1)
}
Ddiff.Monthly <- function(dte1, dte2, dating=NA) {
  stopifnot(inherits(dte1, c("Date", "POSIXt")))
  stopifnot(inherits(dte2, c("Date", "POSIXt")))
  dt1 <- floor_date(dte1, "month")
  dt2 <- floor_date(dte2, "month")
  (year(dt2)-year(dt1))*12+month(dt2)-month(dt1)
}
Ddiff.Weekly <- function(dte1, dte2, dating=NA) {
  stopifnot(inherits(dte1, c("Date", "POSIXt")))
  stopifnot(inherits(dte2, c("Date", "POSIXt")))
  dur <- as.duration(new_interval(floor_date(dte1, "week"), floor_date(dte2, "week")))
  dur@.Data/(3600*24*7)
}
Ddiff.Daily <- function(dte1, dte2, dating=NA) {
  stopifnot(inherits(dte1, c("Date", "POSIXt")))
  stopifnot(inherits(dte2, c("Date", "POSIXt")))
  dur <- as.duration(new_interval(floor_date(dte1, "day"), floor_date(dte2, "day")))
  dur@.Data/(3600*24)
}

#-----------------------------------------------------------------------------
