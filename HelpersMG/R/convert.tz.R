#' convert.tz Convert one Date-Time from one timezone to another
#' @title Convert one Date-Time from one timezone to another
#' @author Marc Girondot
#' @return A POSIXlt or POSIXct date converted
#' @param x The date-time in POSIXlt or POSIXct format
#' @param tz The timezone
#' @description Convert one Date-Time from one timezone to another.\cr
#' Available timezones can be shown using OlsonNames().
#' @seealso Function \code{with_tz()} from \code{lubridate} package does the same. I keep it here only for compatibility with old scripts.
#' @examples
#' d <- as.POSIXlt("2010-01-01 17:34:20", tz="UTC")
#' convert.tz(d, tz="America/Guatemala")
#' @export


convert.tz <- function(x, tz=Sys.timezone()) {
d <- as.POSIXct(format(as.POSIXct(x), tz=tz, usetz=TRUE), tz=tz)
if (any(class(c)=="POSIXlt")) d <- as.POSIXlt(d)
return(d)
}
