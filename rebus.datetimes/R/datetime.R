#' @include imports.R
NULL

#' Date-time regexes
#'
#' Compound regex constants for matching ISO 8601 dates and times.
#' @param x A \code{\link[base]{strptime}}-style date-time format string.
#' @param locale A string specifying a locale.
#' @param io Are you trying to match output or input?  The latter is less
#' strict about leading zeroes and spaces.
#' @return A character vector representing part or all of a regular expression.
#' @note \code{"\%O[dHImMUVwWy]"}, \code{"\%E[cCyYxX]"}, \code{"\%x"}, \code{"\%X"}
#' and \code{"\%+"} are supposed to be locale-dependent upon output,
#' but implementing this in an OS-portable way seems to be much more effort
#' than it's worth.
#' @seealso \code{\link[base]{strptime}} that describes formatting codes,
#' \code{\link{ClassGroups}}, \code{\link[base]{Sys.setlocale}}
#' @examples
#' datetime("%m/%d/%Y")             # match US style dates
#' twelve_or_twentyfour <- rebus.base::or("%H", "%I%p")
#' datetime(twelve_or_twentyfour) # match hours in 24h or 12h format
#'
#' \dontrun{
#' # week days and months can be matched in any locale
#' if(.Platform$OS.type == "windows")
#' {
#'   fr_FR <- "French_France"
#'   ar_QA <- "Arabic_Qatar"
#' } else if(Sys.info()["sysname"] != "Darwin") # mac
#' {
#'   fr_FR <- "fr_FR"
#'   ar_QA <- "ar_QA"
#' } else if(Sys.info()["sysname"] != "Linux")
#' {
#'   fr_FR <- "fr_FR.utf8"
#'   ar_QA <- "ar_QA.utf8"
#' }
#' datetime("%a %A %b %B", fr_FR)
#' datetime("%a %A %b %B", ar_QA)
#'
#' # All letter tokens.  Lots of output.
#' x <- paste0("%", c(letters, LETTERS))
#' stats::setNames(datetime(x), x)
#' }
#'
#' # Individual date-time components
#' DTSEP             # optional selected punctuation or space
#' CENTURY           # exactly two digits
#' YEAR              # one to four digits
#' YEAR2             # exactly two digits
#' YEAR4             # exactly four digits
#' MONTH             # number from 1 to 12, leading zero
#' WEEK_OF_YEAR      # number from 0 to 53, leading zero
#' DAY               # number from 1 to 31, leading zero
#' DAY_SINGLE        # leading space
#' HOUR24            # 24 hour clock, leading zero
#' HOUR12            # 12 hour clock, leading zero
#' HOUR24_SINGLE     # 24 hour clock, leading space
#' HOUR12_SINGLE     # 12 hour clock, leading space
#' MINUTE            # number from 0 to 59, leading zero
#' SECOND            # number from 0 to 61 (leap seconds), leading zero
#' FRACTIONAL_SECOND # a second optional decimal point and up to 6 digits
#' AM_PM             # AM or PM, any case
#' TIMEZONE_OFFSET   # optional plus or minus, then four digits
#' TIMEZONE          # Any value returned by OlsonNames()
#' # ISO 8601 formats
#' ISO_DATE          # %Y-%m-%d
#' ISO_TIME          # %H:%M:%S
#' ISO_DATETIME      # ISO_DATE followed by ISO_TIME, separated by space or "T".
#' # Compound forms, separated by DTSEP
#' YMD
#' YDM
#' MYD
#' MDY
#' DYM
#' DMY
#' HMS
#' HM
#' MS
#'
#' # Some forms have less strict alternatives for input (with an '_IN' suffix).
#' CENTURY_IN
#' MONTH_IN
#' WEEK_OF_YEAR_IN
#' DAY_IN
#' HOUR24_IN
#' HOUR12_IN
#' MINUTE_IN
#' SECOND_IN
#' FRACTIONAL_SECOND_IN
#' ISO_DATE_IN
#' ISO_TIME_IN
#' ISO_DATETIME_IN
#' YMD_IN
#' YDM_IN
#' MYD_IN
#' MDY_IN
#' DYM_IN
#' DMY_IN
#' HMS_IN
#' HM_IN
#' MS_IN
#'
#' dates <- seq(as.Date("2000-01-01"), as.Date("2001-01-01"), "1 day")
#' datetimes <- seq(as.POSIXct(Sys.Date()), as.POSIXct(Sys.Date() + 1), "1 sec")
#' times <- substring(datetimes, 12, 19)
#' stopifnot(
#'   all(grepl(ISO_DATE, dates)),
#'   all(grepl(ISO_TIME, times)),
#'   all(grepl(ISO_DATETIME, datetimes))
#' )
#' non_dates <- c(
#'   "2000-13-01", "2000-01-32", "2000-00-01", "2000-01-00"
#' )
#' non_times <- c(
#'   "24:00:00", "23:60:59", "23:59:62", "23 59 59"
#' )
#' stopifnot(
#'   all(!grepl(ISO_DATE, non_dates)),
#'   all(!grepl(ISO_TIME, non_times))
#' )
#' @name DateTime
NULL

#' Get the days of the week
#'
#' Get the names of the days of the week in a given locale.
#' @param abbreviate A logical value indicating whether or not to return
#' abbreviated names (if they are available).
#' @param locale A string denoting the name of the locale to retrieve names in,
#' or \code{NULL} to use the current locale.
#' @param from Date to use for the first day/month.
#' @return A string containing a regular expression of the names of the days of
#' the week, separated by pipes.  The first day of the week will be the current
#' day.
#' @note See \code{\link[base]{Sys.setlocale}} and
#' \url{http://stackoverflow.com/q/20960821/134830} and
#' \url{http://stackoverflow.com/q/26603564/134830} for how to specify
#' the locale.
#' @examples
#' get_weekdays()
#' get_weekdays(TRUE)
#' get_months()
#' get_months(TRUE)
#' \dontrun{
#' if(.Platform$OS.type == "windows")
#' {
#'   get_weekdays(locale = "French_France")
#'   get_weekdays(TRUE, locale = "French_France")
#'   get_weekdays(locale = "Arabic_Qatar")
#'   get_weekdays(TRUE, locale = "Arabic_Qatar")
#' } else
#' {
#'   get_weekdays(locale = "fr_FR.utf8")
#'   get_weekdays(TRUE, locale = "fr_FR.utf8")
#'   get_weekdays(locale = "ar_QA.utf8")
#'   get_weekdays(TRUE, locale = "ar_QA.utf8")
#' }
#' }
#' @export
get_weekdays <- function(abbreviate = FALSE, locale = NULL, from = Sys.Date())
{
  weekday_names <- if(.Platform$OS.type == "windows")
  {
    get_weekdays_windows(abbreviate, locale, from = from)
  } else
  {
    get_weekdays_posix(abbreviate, locale, from = from)
  }
  or1(weekday_names)
}

get_weekdays_windows <- function(abbreviate = FALSE, locale = NULL, from = Sys.Date())
{
  # TODO: Maybe worth reimplementing devtools::with_locale to clean this up.
  if(!is.null(locale))
  {
    lc_ctype <- Sys.getlocale("LC_CTYPE")
    lc_time <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_CTYPE", lc_ctype))
    on.exit(Sys.setlocale("LC_TIME", lc_time), add = TRUE)
    Sys.setlocale("LC_CTYPE", locale)
    Sys.setlocale("LC_TIME", locale)
  }
  # All locales have 7 days in a week
  days <- seq(from, from + 6, "1 day")
  weekday_names <- weekdays(days, abbreviate)

  # TODO: Is charset <- localeToCharset(Sys.getlocale("LC_TIME"))[1] platform independent?
  current_codepage <- as.character(l10n_info()$codepage)
  iconv(weekday_names, from = current_codepage, to = "utf8")
}

get_weekdays_posix <- function(abbreviate = FALSE, locale = NULL, from = Sys.Date())
{
  if(!is.null(locale))
  {
    lc_time <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_TIME", lc_time), add = TRUE)
    Sys.setlocale("LC_TIME", locale)
  }
  days <- seq(from, from + 6, "1 day")
  weekdays(days, abbreviate)
}

#' @rdname get_weekdays
#' @export
get_months <- function(abbreviate = FALSE, locale = NULL, from = Sys.Date())
{
  month_names <- if(.Platform$OS.type == "windows")
  {
    get_months_windows(abbreviate, locale, from = from)
  } else
  {
    get_months_posix(abbreviate, locale, from = from)
  }
  or1(month_names)
}

get_months_windows <- function(abbreviate = FALSE, locale = NULL, from = Sys.Date())
{
  if(!is.null(locale))
  {
    lc_ctype <- Sys.getlocale("LC_CTYPE")
    lc_time <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_CTYPE", lc_ctype))
    on.exit(Sys.setlocale("LC_TIME", lc_time), add = TRUE)
    Sys.setlocale("LC_CTYPE", locale)
    Sys.setlocale("LC_TIME", locale)
  }
  # All locales have 12 months in a year
  days <- seq(from, by = "1 month", length.out = 12)
  month_names <- months(days, abbreviate)
  current_codepage <- as.character(l10n_info()$codepage)
  iconv(month_names, from = current_codepage, to = "utf8")
}

get_months_posix <- function(abbreviate = FALSE, locale = NULL, from = Sys.Date())
{
  if(!is.null(locale))
  {
    lc_time <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_TIME", lc_time), add = TRUE)
    Sys.setlocale("LC_TIME", locale)
  }
  days <- seq(from, by = "1 month", length.out = 12)
  months(days, abbreviate)
}



#' @rdname DateTime
#' @export
datetime <- function(x, locale = NULL, io = c("output", "input"))
{
  io <- match.arg(io)
  x <- gsub("%E?c", "%a %b %e %H:%M:%S %Y", x)
  x <- gsub("%\\+", "%a %b %e %H:%M:%S %Z %Y", x)
  x <- gsub("%D", "%m/%d/%y", x)
  x <- gsub("%R", "%H:%M", x)
  x <- gsub("%E?x", "%y/%m/%d", x) # This should really be locale-specific on output
  x <- gsub("%[pP]", AM_PM, x)
  x <- gsub("%s", ascii_digit(1, 17), x)
  x <- gsub("%u", WEEKDAY1, x)
  x <- gsub("%O?w", WEEKDAY0, x)
  x <- gsub("%[EO]?y", YEAR2, x)
  x <- gsub("%E?Y", YEAR4, x)
  x <- gsub("%z", TIMEZONE_OFFSET, x)
  x <- gsub("%Z", TIMEZONE, x)
  if(io == "output")
  {
    x <- gsub("%a", get_weekdays(TRUE, locale = locale), x)
    x <- gsub("%A", get_weekdays(locale = locale), x)
    x <- gsub("%[bh]", get_months(TRUE, locale = locale), x)
    x <- gsub("%B", get_months(locale = locale), x)
    x <- gsub("%E?C", CENTURY, x)
    x <- gsub("%O?d", DAY, x)
    x <- gsub("%e", DAY_SINGLE, x)
    x <- gsub("%F", ISO_DATE, x)
    x <- gsub("%O?H", HOUR24, x)
    x <- gsub("%O?[Ir]", HOUR12, x)
    x <- gsub("%j", DAY_OF_YEAR, x)
    x <- gsub("%k", HOUR24_SINGLE, x)
    x <- gsub("%l", HOUR12_SINGLE, x)
    x <- gsub("%O?m", MONTH, x)
    x <- gsub("%O?M", MINUTE, x)
    x <- gsub("%n", "\\n", x)
    x <- gsub("%OS", FRACTIONAL_SECOND, x)
    x <- gsub("%S", SECOND, x)
    x <- gsub("%t", "\\t", x)
    x <- gsub("%(?:T|E?X)", ISO_TIME, x) # %X should really be locale-specific on output
    x <- gsub("%O?[gGUVW]", WEEK_OF_YEAR, x)
  } else # io == "input"
  {
    x <- gsub(
      "%aA",
      or(
        get_weekdays(TRUE, locale = locale),
        get_weekdays(locale = locale)
      ),
      x
    )
    x <- gsub(
      "%[bBh]",
      or(
        get_months(TRUE, locale = locale),
        get_months(locale = locale)
      ),
      x
    )
    x <- gsub("%E?C", CENTURY_IN, x)
    x <- gsub("%(?:O?d|e)", DAY_IN, x)
    x <- gsub("%F", ISO_DATE_IN, x)
    x <- gsub("%(?:O?H|k)", HOUR24_IN, x)
    x <- gsub("%(?:O?[Ir]|l)", HOUR12_IN, x)
    x <- gsub("%j", DAY_OF_YEAR_IN, x)
    x <- gsub("%O?m", MONTH_IN, x)
    x <- gsub("%O?M", MINUTE_IN, x)
    x <- gsub("%[nt]", space(0, Inf), x)
    x <- gsub("%OS", FRACTIONAL_SECOND_IN, x)
    x <- gsub("%S", SECOND_IN, x)
    x <- gsub("%(?:T|E?X)", ISO_TIME_IN, x)
    x <- gsub("%O?[gGUVW]", WEEK_OF_YEAR_IN, x)
    x <- case_insensitive(x)
  }
  group(x)
}
