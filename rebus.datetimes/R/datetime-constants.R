#' @include imports.R
NULL

#' @rdname DateTime
#' @export
OPT_LEADING_0 <- "[0 ]?"

#' @rdname DateTime
#' @importFrom rebus.base optional
#' @importFrom rebus.base char_class
#' @export
DTSEP <- optional(char_class("-/.:,\\ "))

#' @rdname DateTime
#' @export
CENTURY <- ascii_digit(2)

#' @rdname DateTime
#' @export
CENTURY_IN <- group(
  OPT_LEADING_0 %R% ascii_digit() %|%
    ascii_digit(2)
)

#' @rdname DateTime
#' @export
YEAR <- ascii_digit(1, 4)

#' @rdname DateTime
#' @export
YEAR2 <- ascii_digit(2)

#' @rdname DateTime
#' @export
YEAR4 <- ascii_digit(4)

#' @rdname DateTime
#' @export
MONTH <- group(
  "0" %R% char_range(1, 9) %|%
    "1" %R% char_range(0, 2)
)

#' @rdname DateTime
#' @export
MONTH_IN <- group(
  OPT_LEADING_0 %R% char_range(1, 9) %|%
    "1" %R% char_range(0, 2)
)

#' @rdname DateTime
#' @export
WEEK_OF_YEAR <- group(
  char_range(0, 4) %R% ascii_digit() %|%
    "5" %R% char_range(0, 3)
)

#' @rdname DateTime
#' @export
WEEK_OF_YEAR_IN <- group(
  OPT_LEADING_0 %R% ascii_digit() %|%
  char_range(1, 4) %R% ascii_digit() %|%
    "5" %R% char_range(0, 3)
)

#' @rdname DateTime
#' @export
DAY <- group(
  "0" %R% char_range(1, 9) %|%
    char_class("12") %R% ascii_digit() %|%
    "3" %R% char_class("01")
)

#' @rdname DateTime
#' @export
DAY_IN <- group(
  OPT_LEADING_0 %R% char_range(1, 9) %|%
    char_class("12") %R% ascii_digit() %|%
    "3" %R% char_class("01")
)

#' @rdname DateTime
#' @export
DAY_SINGLE <- group(
  " " %R% char_range(1, 9) %|%
    char_class("12") %R% ascii_digit() %|%
    "3" %R% char_class("01")
)

#' @rdname DateTime
#' @export
DAY_OF_YEAR <- group(
  "00" %R% char_range(1, 9) %|%
    "0" %R% char_range(1, 9) %R% ascii_digit() %|%
    char_class("12") %R% ascii_digit(2) %|%
    "3" %R% char_range(0, 5) %R% ascii_digit() %|%
    "36" %R% char_range(0, 6)
)

#' @rdname DateTime
#' @export
DAY_OF_YEAR_IN <- group(
  repeated(char_class("0 "), 0, 2) %R% char_range(1, 9) %|%
    OPT_LEADING_0 %R% char_range(1, 9) %R% ascii_digit() %|%
    char_class("12") %R% ascii_digit(2) %|%
    "3" %R% char_range(0, 5) %R% ascii_digit() %|%
    "36" %R% char_range(0, 6)
)

#' @rdname DateTime
#' @export
WEEKDAY1 <- char_range(1, 7)

#' @rdname DateTime
#' @export
WEEKDAY0 <- char_range(0, 6)

#' @rdname DateTime
#' @export
HOUR24 <- group(
  char_class("01") %R% ascii_digit() %|% "2" %R% char_range("0", "3")
)

#' @rdname DateTime
#' @export
HOUR24_SINGLE <- group(
  char_class(" 1") %R% ascii_digit() %|% "2" %R% char_range("0", "3")
)

#' @rdname DateTime
#' @export
HOUR24_IN <- group(
  optional(char_class(" 01")) %R% ascii_digit() %|%
    "2" %R% char_range("0", "4") # 24:00:00 allowed on input
)

#' @rdname DateTime
#' @export
HOUR12 <- group(
  "0" %R% char_range("1", "9")  %|% "1" %R% char_range("0", "2")
)

#' @rdname DateTime
#' @export
HOUR12_SINGLE <- group(
  " " %R% char_range("1", "9") %|% "1" %R% char_range("0", "2")
)

#' @rdname DateTime
#' @export
HOUR12_IN <- group(
  OPT_LEADING_0 %R% char_range("1", "9") %|% "1" %R% char_range("0", "2")
)

#' @rdname DateTime
#' @export
MINUTE <- char_range(0, 5) %R% ascii_digit()

#' @rdname DateTime
#' @export
MINUTE_IN <- optional(char_class(" 0-5")) %R% ascii_digit()

#' @rdname DateTime
#' @export
SECOND <- group(
  char_range(0, 5) %R% ascii_digit() %|%
  "6" %R% char_class("01") # leap seconds
)

#' @rdname DateTime
#' @export
SECOND_IN <- group(
  optional(char_class(" 0-5")) %R% ascii_digit() %|%
    "6" %R% char_class("01") # leap seconds
)

#' @rdname DateTime
#' @export
FRACTIONAL_SECOND <- SECOND %R%
  optional(group(char_class(".,") %R% ascii_digit(1, 6)))

#' @rdname DateTime
#' @export
FRACTIONAL_SECOND_IN <- SECOND_IN %R%
  optional(group(char_class(".,") %R% ascii_digit(1, 6)))

#' @rdname DateTime
#' @export
AM_PM <- group("am" %|% "AM" %|% "pm" %|% "PM")

#' @rdname DateTime
#' @export
TIMEZONE_OFFSET <- optional(char_class("-+")) %R% ascii_digit(4)

#' @rdname DateTime
#' @export
TIMEZONE <- or1(escape_special(OlsonNames()))


#' @rdname DateTime
#' @export
ISO_DATE <- YEAR4 %R% "-" %R% MONTH %R% "-" %R% DAY

#' @rdname DateTime
#' @export
ISO_DATE_IN <- YEAR4 %R% "-" %R% MONTH_IN %R% "-" %R% DAY_IN

#' @rdname DateTime
#' @export
ISO_TIME <- HOUR24 %R%  ":" %R% MINUTE %R%  ":" %R% SECOND

#' @rdname DateTime
#' @export
ISO_TIME_IN <- HOUR24_IN %R%  ":" %R% MINUTE_IN %R%  ":" %R% SECOND_IN

#' @rdname DateTime
#' @export
ISO_DATETIME <- ISO_DATE %R% char_class(" T") %R% ISO_TIME

#' @rdname DateTime
#' @export
ISO_DATETIME_IN <- ISO_DATE_IN %R% char_class(" T") %R% ISO_TIME_IN


#' @rdname DateTime
#' @export
YMD <- YEAR %R% DTSEP %R% MONTH %R% DTSEP %R% DAY

#' @rdname DateTime
#' @export
YMD_IN <- YEAR %R% DTSEP %R% MONTH_IN %R% DTSEP %R% DAY_IN

#' @rdname DateTime
#' @export
YDM <- YEAR %R% DTSEP %R% DAY %R% DTSEP %R% MONTH

#' @rdname DateTime
#' @export
YDM_IN <- YEAR %R% DTSEP %R% DAY_IN %R% DTSEP %R% MONTH_IN

#' @rdname DateTime
#' @export
MYD <- MONTH %R% DTSEP %R% YEAR %R% DTSEP %R% DAY

#' @rdname DateTime
#' @export
MYD_IN <- MONTH_IN %R% DTSEP %R% YEAR %R% DTSEP %R% DAY_IN

#' @rdname DateTime
#' @export
MDY <- MONTH %R% DTSEP %R% DAY %R% DTSEP %R% YEAR

#' @rdname DateTime
#' @export
MDY_IN <- MONTH_IN %R% DTSEP %R% DAY_IN %R% DTSEP %R% YEAR

#' @rdname DateTime
#' @export
DYM <- DAY %R% DTSEP %R% YEAR %R% DTSEP %R% MONTH

#' @rdname DateTime
#' @export
DYM_IN <- DAY_IN %R% DTSEP %R% YEAR %R% DTSEP %R% MONTH_IN

#' @rdname DateTime
#' @export
DMY <- DAY %R% DTSEP %R% MONTH %R% DTSEP %R% YEAR

#' @rdname DateTime
#' @export
DMY_IN <- DAY_IN %R% DTSEP %R% MONTH_IN %R% DTSEP %R% YEAR

#' @rdname DateTime
#' @export
HMS <- HOUR24 %R% DTSEP %R% MINUTE %R% DTSEP %R% SECOND

#' @rdname DateTime
#' @export
HMS_IN <- HOUR24_IN %R% DTSEP %R% MINUTE_IN %R% DTSEP %R% SECOND_IN

#' @rdname DateTime
#' @export
HM <- HOUR24 %R% DTSEP %R% MINUTE

#' @rdname DateTime
#' @export
HM_IN <- HOUR24_IN %R% DTSEP %R% MINUTE_IN

#' @rdname DateTime
#' @export
MS <- MINUTE %R% DTSEP %R% SECOND

#' @rdname DateTime
#' @export
MS_IN <- MINUTE_IN %R% DTSEP %R% SECOND_IN
