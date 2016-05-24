
#-----------------------------------------------------------------------------
# Weekdays

Mondays <- new_Dating(c("Weekdays", "Mondays"))
Tuesdays <- new_Dating(c("Weekdays", "Tuesdays"))
Wednesdays <- new_Dating(c("Weekdays", "Wednesdays"))
Thursdays <- new_Dating(c("Weekdays", "Thursdays"))
Fridays <- new_Dating(c("Weekdays", "Fridays"))
Saturdays <- new_Dating(c("Weekdays", "Saturdays"))
Sundays <- new_Dating(c("Weekdays", "Sundays"))

.DDweekday <- function(x) UseMethod(".DDweekday")

.DDweekday.Sundays <- function(x) { 0 }
.DDweekday.Mondays <- function(x) { 1 }
.DDweekday.Tuesdays <- function(x) { 2 }
.DDweekday.Wednesdays <- function(x) { 3 }
.DDweekday.Thursdays <- function(x) { 4 }
.DDweekday.Fridays <- function(x) { 5 }
.DDweekday.Saturdays <- function(x) { 6 }

Dbelong.Weekdays <- function(dte, dating) {
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  stopifnot(inherits(dating, "Weekdays"))
  wd <- .DDweekday(dating)
  if(wd>0) Dbelong(Dsucc(dte,Daily,-wd), Weekly)
  else Dbelong(dte, Weekly)
}

Dseq.Weekdays <- function(from, to, dating, len) {
  stopifnot(inherits(dating, "Weekdays"))
  stopifnot(inherits(from, c("Date", "POSIXt")))
  from <- as.Date(from)
  wd <- .DDweekday(dating)
  if(missing(to)) {
    if(wd>0) Dsucc(Dseq(Dsucc(from,Daily,-wd), , Weekly, len), Daily, wd)
    else Dseq(from, , Weekly, len)
  } else {
    stopifnot(inherits(to, c("Date", "POSIXt")))
    to <- as.Date(to)
    stopifnot(from <= to)
    if(wd>0) Dsucc(Dseq(Dsucc(from,Daily,-wd), Dsucc(to,Daily,-wd), Weekly), Daily, wd)
    else Dseq(from, to, Weekly)
  }
}

Dsucc.Weekdays <- function(dte, dating, num=1) {
  stopifnot(inherits(dating, "Weekdays"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  wd <- .DDweekday(dating)
  if(wd>0) Dsucc(Dsucc(Dsucc(dte, Daily, -wd), Weekly, num), Daily, wd)
  else Dsucc(dte, Weekly, num)
}

Dfloor.Weekdays <- function(dte, dating) {
  stopifnot(inherits(dating, "Weekdays"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  wd <- .DDweekday(dating)
  if(wd>0) Dsucc(Dfloor(Dsucc(dte, Daily, -wd), Weekly), Daily, wd)
  else Dfloor(dte, Weekly)
}

Dceiling.Weekdays <- function(dte, dating) {
  stopifnot(inherits(dating, "Weekdays"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  wd <- .DDweekday(dating)
  if(wd>0) Dsucc(Dceiling(Dsucc(dte, Daily, -wd), Weekly), Daily, wd)
  else Dceiling(dte, Weekly)
}

Dround.Weekdays <- function(dte, dating) {
  stopifnot(inherits(dating, "Weekdays"))
  stopifnot(inherits(dte, c("Date", "POSIXt")))
  dte <- as.Date(dte)
  wd <- .DDweekday(dating)
  if(wd>0) Dsucc(Dround(Dsucc(dte, Daily, -wd), Weekly), Daily, wd)
  else Dround(dte, Weekly)
}

Ddiff.Weekdays <- function(dte1, dte2, dating) {
  stopifnot(inherits(dating, "Weekdays"))
  stopifnot(inherits(dte1, c("Date", "POSIXt")))
  stopifnot(inherits(dte2, c("Date", "POSIXt")))
  dte1 <- as.Date(dte1)
  dte2 <- as.Date(dte2)
  wd <- .DDweekday(dating)
  if(wd>0) Ddiff(Dsucc(dte1, Daily, -wd), Dsucc(dte2, Daily, -wd), Weekly)
  else Ddiff(dte1, dte2, Weekly)
}

#-----------------------------------------------------------------------------
