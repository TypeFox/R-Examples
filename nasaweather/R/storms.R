#' Storm tracks data
#'
#' Tropical cyclone tracks through the Atlantic Ocean, Caribbean Sea and Gulf
#' of Mexico from 1995 to 2005. Only "named" storms, those which reached
#' tropical storm status or stronger, are included.
#'
#' The data originated from the National Hurricane Center's archive of
#' Tropical Cyclone Reports (\url{http://www.nhc.noaa.gov/}). This dataset
#' was hand-scraped from best track tables in the individual tropical cyclone
#' reports (PDF, HTML and Microsoft Word) by Jon Hobbs.
#'
#' The Tropical Cyclone Reports had a variety of storm type designations and
#' there appeared to be no consistent naming convention for cyclones that were
#' not hurricanes, tropical depressions, or tropical storms. Many of these
#' designations have been combined into the "Extratropical" category in this
#' dataset.
#'
#' @author Jon Hobbs
#' @format
#' \describe{
#' \item{name}{Storm Name}
#' \item{year,month,day}{Date of report}
#' \item{hour}{Hour of report (0, 6, 12 or 18, in UTC)}
#' \item{lat,long}{Location of storm center}
#' \item{pressure}{Air pressure at the storm's center (in millibars)}
#' \item{wind}{wtorm's maximum sustained wind speed (in knots)}
#' \item{type}{Storm classification (Tropical Depression, Tropical Storm,
#'   Hurricane, or Extratropical)}
#' \item{seasday}{Day of the hurricane season (days since June 1)}
#' }
"storms"
