### getDateTimePart.R
###------------------------------------------------------------------------
### What: Extract date and time parts from various date and time classes
### $Id$
### Time-stamp: <2008-12-30 22:42:58 ggorjan>
###------------------------------------------------------------------------

### {{{ getYear
###------------------------------------------------------------------------

getYear <- function(x, format, ...)
  UseMethod("getYear")

getYear.default <- function(x, format, ...)
  stop("'getYear' can only be used on objects of a date/time class")

getYear.Date <-
getYear.POSIXct <-
getYear.POSIXlt <- function(x, format="%Y", ...)
  format(x=x, format=format, ...)

### }}}
### {{{ getMonth
###------------------------------------------------------------------------

getMonth <- function(x, format, ...)
  UseMethod("getMonth")

getMonth.default <- function(x, format, ...)
  stop("'getMonth' can only be used on objects of a date/time class")

getMonth.Date <-
getMonth.POSIXct <-
getMonth.POSIXlt <- function(x, format="%m", ...)
  format(x=x, format=format)

### }}}
### {{{ getDay
###------------------------------------------------------------------------

getDay <- function(x, format, ...)
  UseMethod("getDay")

getDay.default <- function(x, format, ...)
  stop("'getDay' can only be used on objects of a date/time class")

getDay.Date <-
getDay.POSIXct <-
getDay.POSIXlt <- function(x, format="%d", ...)
  format(x=x, format=format)

### }}}
### {{{ getHour
###------------------------------------------------------------------------

getHour <- function(x, format, ...)
  UseMethod("getHour")

getHour.default <- function(x, format, ...)
  stop("'getHour' can only be used on objects of a date/time class")

### }}}
### {{{ getMin
###------------------------------------------------------------------------

getMin <- function(x, format, ...)
  UseMethod("getMin")

getMin.default <- function(x, format, ...)
  stop("'getMin' can only be used on objects of a date/time class")

### }}}
### {{{ getSec
###------------------------------------------------------------------------

getSec <- function(x, format, ...)
  UseMethod("getSec")

getSec.default <- function(x, format, ...)
  stop("'getSec' can only be used on objects of a date/time class")

### }}}
### {{{ Dear Emacs
## Local variables:
## folded-file: t
## End:
### }}}

###------------------------------------------------------------------------
### getDateTimePart.R ends here